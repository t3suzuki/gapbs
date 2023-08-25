#include <iostream>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "generator");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  Writer w(g);
  w.WriteGraph("new.sg", true);
  return 0;
}
