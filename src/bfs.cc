// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <iostream>
#include <vector>

#include <signal.h>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"
#include "dax.h"

/*
GAP Benchmark Suite
Kernel: Breadth-First Search (BFS)
Author: Scott Beamer

Will return parent array for a BFS traversal from a source vertex

This BFS implementation makes use of the Direction-Optimizing approach [1].
It uses the alpha and beta parameters to determine whether to switch search
directions. For representing the frontier, it uses a SlidingQueue for the
top-down approach and a Bitmap for the bottom-up approach. To reduce
false-sharing for the top-down approach, thread-local QueueBuffer's are used.

To save time computing the number of edges exiting the frontier, this
implementation precomputes the degrees in bulk at the beginning by storing
them in parent array as negative numbers. Thus the encoding of parent is:
  parent[x] < 0 implies x is unvisited and parent[x] = -out_degree(x)
  parent[x] >= 0 implies x been visited

[1] Scott Beamer, Krste Asanović, and David Patterson. "Direction-Optimizing
    Breadth-First Search." International Conference on High Performance
    Computing, Networking, Storage and Analysis (SC), Salt Lake City, Utah,
    November 2012.
*/

//#define PTH (1)

#define STRIDE (128)
#define N_TH (8)
#define ULT_N_TH (8*N_TH)

#if PTH
#include <pth.h>
#endif

#if ABT
#include <abt.h>
static ABT_xstream abt_xstreams[N_TH];
static ABT_thread abt_threads[ULT_N_TH];
static ABT_pool abt_pools[N_TH];
#endif

using namespace std;

int64_t BUStep(const Graph &g, pvector<NodeID> &parent, Bitmap &front,
               Bitmap &next) {
  int64_t awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    if (parent[u] < 0) {
      for (NodeID v : g.in_neigh(u)) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}


int64_t TDStep(const Graph &g, pvector<NodeID> &parent,
               SlidingQueue<NodeID> &queue) {
  int64_t scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
#pragma omp for reduction(+ : scout_count) nowait schedule(dynamic, 64)
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (NodeID v : g.out_neigh(u)) {
        NodeID curr_val = parent[v];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.push_back(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}




bool
__get_next(size_t &qp,
	   SlidingQueue<NodeID> &queue,
	   SlidingQueue<NodeID>::iterator &q_iter,
	   SlidingQueue<NodeID>::iterator &q_iter_end,
	   int stride)
{
  size_t old_qp = fetch_and_add(qp, stride);
  size_t new_qp = old_qp + stride;
  if (old_qp < queue.size()) {
    q_iter = queue.begin() + old_qp;
    q_iter_end = (new_qp >= queue.size()) ? queue.end() : queue.begin() + new_qp;
    return true;
  } else {
    q_iter = queue.end();
    q_iter_end = queue.end();
    return false;
  }
}

inline bool get_next(size_t &qp,
		     SlidingQueue<NodeID> &queue,
		     SlidingQueue<NodeID>::iterator &q_iter,
		     SlidingQueue<NodeID>::iterator &q_iter_end)
{
  return __get_next(qp, queue, q_iter, q_iter_end, STRIDE);
}

struct arg_t {
  int id;
  const Graph *g;
  pvector<NodeID> *parent;
  SlidingQueue<NodeID> *queue;
  size_t *qp;
  int64_t *scout_count;
  SlidingQueue<NodeID>::iterator q_iter;
  SlidingQueue<NodeID>::iterator q_iter_end;
#if PTH
  pth_t pthid;
#endif
};

#define N_ENTRY (1024/ULT_N_TH)
uint64_t mycache[ULT_N_TH][N_ENTRY];

inline char *align(void *p, int sz)
 {
  return (char *)(((uint64_t)p / sz) * sz);
}

void
thread_func(arg_t *arg)
{
  SlidingQueue<NodeID>::iterator q_iter;
  SlidingQueue<NodeID>::iterator q_iter_end;
  QueueBuffer<NodeID> lqueue(*arg->queue);
  int64_t local_scout_count = 0;
  while (get_next(*arg->qp, *arg->queue, q_iter, q_iter_end)) {
    while (q_iter < q_iter_end) {
      NodeID u = *q_iter++;
#if ABT
      /*
      void *p = arg->g->out_top(u);
      __builtin_prefetch(p);
      ABT_thread_yield();
      */
#define ALIGN (64)
      if (0) {
	char *p = align(arg->g->out_top(u), ALIGN);
	int tid = arg->id;
	int index = ((uint64_t)p / ALIGN) % N_ENTRY;
	if (mycache[tid][index] != (uint64_t)p / ALIGN) {
	  for (int i=0; i<ALIGN/64; i++) {
	    char *np = p + 64 * i;
	    __builtin_prefetch(np);
	  }
	  mycache[tid][index] = (uint64_t)p / ALIGN;
	  ABT_thread_yield();
	}
      } else {
	arg->g->out_prefetch2(u);
	ABT_thread_yield();
      }
#endif

      for (NodeID v : arg->g->out_neigh(u)) {
	__builtin_prefetch(&(*arg->parent)[v]);
	ABT_thread_yield();
	NodeID curr_val = (*arg->parent)[v];
	if (curr_val < 0) {
	  if (compare_and_swap((*arg->parent)[v], curr_val, u)) {
	    lqueue.push_back(v);
	    local_scout_count += -curr_val;
	  }
	}
      }
    }
  }
  fetch_and_add(*arg->scout_count, local_scout_count);
  lqueue.flush();
}

int64_t TDStep_pthread(const Graph &g, pvector<NodeID> &parent,
		       SlidingQueue<NodeID> &queue) {
  int64_t scout_count = 0;
  size_t qp = 0;
  pthread_t pth[N_TH];
  arg_t arg[N_TH];
  for (int i=0; i<N_TH; i++) {
    arg[i].id = i;
    arg[i].g = &g;
    arg[i].parent = &parent;
    arg[i].queue = &queue;
    arg[i].qp = &qp;
    arg[i].scout_count = &scout_count;
    //pthread_setname_np(pth[i], "worker");
    pthread_create(&pth[i], NULL, (void * (*)(void *))thread_func, &arg[i]);
  }
  for (int i=0; i<N_TH; i++) {
    pthread_join(pth[i], NULL);
  }
  return scout_count;
}


#if ABT
void
my_join(void *abt_th)
{
  ABT_thread_free((ABT_thread *)abt_th);
}

int64_t TDStep_abt(const Graph &g, pvector<NodeID> &parent,
		       SlidingQueue<NodeID> &queue) {
  int64_t scout_count = 0;
  size_t qp = 0;
  arg_t arg[ULT_N_TH];
  int i;
  for (i=0; i<ULT_N_TH; i++) {
    arg[i].id = i;
    arg[i].g = &g;
    arg[i].parent = &parent;
    arg[i].queue = &queue;
    arg[i].qp = &qp;
    arg[i].scout_count = &scout_count;
    ABT_thread_create(abt_pools[i % N_TH], (void (*)(void *))thread_func, &arg[i], ABT_THREAD_ATTR_NULL, &abt_threads[i]);
  }
  for (i=0; i<ULT_N_TH; i++) {
    ABT_thread_join(abt_threads[i]);
  }
  return scout_count;
}
#endif

#if PTH
int64_t TDStep_pth(const Graph &g, pvector<NodeID> &parent,
		   SlidingQueue<NodeID> &queue) {
  int64_t scout_count = 0;
  size_t qp = 0;
  arg_t arg[ULT_N_TH];
  int i;
  pth_attr_t attr = pth_attr_new();
  for (i=0; i<ULT_N_TH; i++) {
    arg[i].id = i;
    arg[i].g = &g;
    arg[i].parent = &parent;
    arg[i].queue = &queue;
    arg[i].qp = &qp;
    arg[i].scout_count = &scout_count;
    arg[i].pthid = pth_spawn(attr, (void *(*)(void *))thread_func, &arg[i]);
  }
  for (i=0; i<ULT_N_TH; i++) {
    pth_join(arg[i].pthid, NULL);
  }
  return scout_count;
}
#endif

void
thread_func2(arg_t *arg)
{
  QueueBuffer<NodeID> lqueue(*(arg->queue));
  const Graph *pg = arg->g;
  pvector<NodeID> *pparent = arg->parent;
  int64_t local_scout_count = 0;
  SlidingQueue<NodeID>::iterator iter_start = arg->q_iter;
  SlidingQueue<NodeID>::iterator iter_end = arg->q_iter_end;
  for (auto q_iter = iter_start; q_iter < iter_end; q_iter++) {
    NodeID u = *q_iter;
    for (NodeID v : pg->out_neigh(u)) {
      NodeID curr_val = (*pparent)[v];
      if (curr_val < 0) {
        if (compare_and_swap((*pparent)[v], curr_val, u)) {
          lqueue.push_back(v);
          local_scout_count += -curr_val;
        }
      }
    }
  }
  fetch_and_add(*arg->scout_count, local_scout_count);
  lqueue.flush();
}

int64_t TDStep_pthread2(const Graph &g, pvector<NodeID> &parent,
			SlidingQueue<NodeID> &queue) {
  int64_t scout_count = 0;
  pthread_t pth[N_TH];
  arg_t arg[N_TH];

  int stride = queue.size() / N_TH;
  for (int tid=0; tid<N_TH; tid++) {
    arg[tid].q_iter = queue.begin() + stride * tid;
    arg[tid].q_iter_end = queue.begin() + stride * (tid + 1);
  }
  arg[N_TH-1].q_iter_end = queue.end();

  for (int i=0; i<N_TH; i++) {
    arg[i].g = &g;
    arg[i].parent = &parent;
    arg[i].queue = &queue;
    arg[i].scout_count = &scout_count;
    pthread_create(&pth[i], NULL, (void * (*)(void *))thread_func2, &arg[i]);
  }
  for (int i=0; i<N_TH; i++) {
    pthread_join(pth[i], NULL);
  }  
  return scout_count;
}


void QueueToBitmap(const SlidingQueue<NodeID> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

void BitmapToQueue(const Graph &g, const Bitmap &bm,
                   SlidingQueue<NodeID> &queue) {
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for nowait
    for (NodeID n=0; n < g.num_nodes(); n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  queue.slide_window();
}

pvector<NodeID> InitParent(const Graph &g) {
  pvector<NodeID> parent(g.num_nodes());
  #pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    parent[n] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}

pvector<NodeID> __DOBFS(const Graph &g, NodeID source, int alpha = 15,
                      int beta = 18) {
  PrintStep("Source", static_cast<int64_t>(source));
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g);
  t.Stop();
  PrintStep("i", t.Seconds());
  parent[source] = source;
  SlidingQueue<NodeID> queue(g.num_nodes());
  queue.push_back(source);
  queue.slide_window();
  Bitmap curr(g.num_nodes());
  curr.reset();
  Bitmap front(g.num_nodes());
  front.reset();
  int64_t edges_to_check = g.num_edges_directed();
  int64_t scout_count = g.out_degree(source);
  while (!queue.empty()) {
    //if (scout_count > edges_to_check / alpha) {
    if (0) {
      int64_t awake_count, old_awake_count;
      TIME_OP(t, QueueToBitmap(queue, front));
      PrintStep("e", t.Seconds());
      awake_count = queue.size();
      queue.slide_window();
      do {
        t.Start();
        old_awake_count = awake_count;
        awake_count = BUStep(g, parent, front, curr);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while ((awake_count >= old_awake_count) ||
               (awake_count > g.num_nodes() / beta));
      TIME_OP(t, BitmapToQueue(g, front, queue));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
#if ABT
      scout_count = TDStep_abt(g, parent, queue);
#elif PTH
      scout_count = TDStep_pth(g, parent, queue);
#else
#if defined _OPENMP
      scout_count = TDStep(g, parent, queue);
#else
      scout_count = TDStep_pthread(g, parent, queue);
#endif
#endif
      queue.slide_window();
      t.Stop();
      PrintStep("td", t.Seconds(), queue.size());
    }
  }
  #pragma omp parallel for
  for (NodeID n = 0; n < g.num_nodes(); n++)
    if (parent[n] < -1)
      parent[n] = -1;
  return parent;
}


pvector<NodeID> DOBFS(const Graph &g, NodeID source, int alpha = 15,
                      int beta = 18) {
#if 0
  int pid= getpid();
  int cpid = fork();
  if( cpid == 0) {
    // child process .  Run your perf stat
    char buf[256];
    //sprintf(buf, "sudo perf stat -ddd -C 0-9 -p %d > stat.log 2>&1", pid);
    //sprintf(buf, "sudo perf stat -ddd -C 0-9 -p %d", pid);
    //sprintf(buf, "sudo perf stat -e mem_load_retired.local_pmm -e LLC-loads -e LLC-load-misses -e L1-dcache-loads -e L1-dcache-load-misses -e cycles -e instructions -e dTLB-loads -e dTLB-load-misses -C 16-31 -p %d", pid);
    sprintf(buf, "sudo perf record -ag");
    execl("/bin/sh", "sh", "-c", buf, NULL);
    while (1);
  } else {
    setpgid(cpid, 0);
    sleep(1);
    pvector<NodeID> ret = __DOBFS(g, source, alpha, beta);
    kill(-cpid, SIGINT);
    sleep(1);
    return ret;
  }
#else
  return __DOBFS(g, source, alpha, beta);
#endif
}


void PrintBFSStats(const Graph &g, const pvector<NodeID> &bfs_tree) {
  int64_t tree_size = 0;
  int64_t n_edges = 0;
  for (NodeID n : g.vertices()) {
    if (bfs_tree[n] >= 0) {
      n_edges += g.out_degree(n);
      tree_size++;
    }
  }
  cout << "BFS Tree has " << tree_size << " nodes and ";
  cout << n_edges << " edges" << endl;
}


// BFS verifier does a serial BFS from same source and asserts:
// - parent[source] = source
// - parent[v] = u  =>  depth[v] = depth[u] + 1 (except for source)
// - parent[v] = u  => there is edge from u to v
// - all vertices reachable from source have a parent
bool BFSVerifier(const Graph &g, NodeID source,
                 const pvector<NodeID> &parent) {
  pvector<int> depth(g.num_nodes(), -1);
  depth[source] = 0;
  vector<NodeID> to_visit;
  to_visit.reserve(g.num_nodes());
  to_visit.push_back(source);
  for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
    NodeID u = *it;
    for (NodeID v : g.out_neigh(u)) {
      if (depth[v] == -1) {
        depth[v] = depth[u] + 1;
        to_visit.push_back(v);
      }
    }
  }
  for (NodeID u : g.vertices()) {
    if ((depth[u] != -1) && (parent[u] != -1)) {
      if (u == source) {
        if (!((parent[u] == u) && (depth[u] == 0))) {
          cout << "Source wrong" << endl;
          return false;
        }
        continue;
      }
      bool parent_found = false;
      for (NodeID v : g.in_neigh(u)) {
        if (v == parent[u]) {
          if (depth[v] != depth[u] - 1) {
            cout << "Wrong depths for " << u << " & " << v << endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        cout << "Couldn't find edge from " << parent[u] << " to " << u << endl;
        return false;
      }
    } else if (depth[u] != parent[u]) {
      cout << "Reachability mismatch" << endl;
      return false;
    }
  }
  return true;
}


int main(int argc, char* argv[]) {
#if DAX
  dax_init();
#endif  
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();

#if ABT
  int i;
  printf("argobots!\n");
  ABT_init(0, NULL);
  ABT_xstream_self(&abt_xstreams[0]);
  for (i=1; i<N_TH; i++) {
    ABT_xstream_create(ABT_SCHED_NULL, &abt_xstreams[i]);
  }
  for (i=0; i<N_TH; i++) {
    ABT_xstream_set_cpubind(abt_xstreams[i], i);
    ABT_xstream_get_main_pools(abt_xstreams[i], 1, &abt_pools[i]);
  }
#endif
#if PTH
  printf("GNU Pth\n");
  pth_init();
#endif

  SourcePicker<Graph> sp(g, cli.start_vertex());
  auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext()); };
  SourcePicker<Graph> vsp(g, cli.start_vertex());
  auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
    return BFSVerifier(g, vsp.PickNext(), parent);
  };
  BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);
  return 0;
}
