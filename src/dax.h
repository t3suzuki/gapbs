#pragma once

#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "memkind.h"
struct memkind *pmem_kind = NULL;

static char *dax_base = NULL;

void
dax_init()
{
#define MMAP_SIZE ((size_t)1024*1024*1024*16)
  char path[] = "/dev/dax0.0";
  int fd = open(path, O_RDWR);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  dax_base = (char *)mmap(NULL, MMAP_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED|MAP_POPULATE, fd, 0);
  int err = memkind_create_fixed(dax_base, MMAP_SIZE, &pmem_kind);
  if (err) {
    exit(-1);
  }

  printf("DAX init done\n");
}

void *
dax_malloc(size_t sz)
{
  return memkind_malloc(pmem_kind, sz);
}
