
#include <stdlib.h>
#include <stdio.h>
#include <stdalign.h>
#include <stdint.h>

#include <time.h>
#include <assert.h>
#include <string.h>

#include <x86intrin.h>

#include <emmintrin.h>

#define TEST_N 3
#define ITER 5000000


static uint64_t now_ms(void) {
  return (uint64_t)(clock() * 1000 / CLOCKS_PER_SEC);
}

static void print_bits_u32(unsigned int s) {
  printf("%X - ", s);  
  for (int b=31;b>=0;b--) {
    printf("%d", s & (1<<b) ? 1:0);
    if (b%8==0)
      fputc(' ',stdout); 
  }
  fputc('\n',stdout); 
}


static void print_bits_u8(unsigned char s) {
  printf("%X - ", s);  
  for (int b=7;b>=0;b--) {
    printf("%d", s & (1<<b) ? 1:0);
    if (b%8==0)
      fputc(' ',stdout); 
  }
  fputc('\n',stdout); 
}


unsigned char minbits_table[256] = {
0,
1,2,2,3,
3,3,3,4,
4,4,4,4,
4,4,4,5,
5,5,5,5,
5,5,5,5,
5,5,5,5,
5,5,5,6,
6,6,6,6,
6,6,6,6,
6,6,6,6,
6,6,6,6,
6,6,6,6,
6,6,6,6,
6,6,6,6,
6,6,6,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,7,
7,7,7,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,8,
8,8,8,};


int sizeofints_og(int num_of_ints, unsigned int sizes[])
{
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;

    /* potential vectorise here */
    for (i=0; i < num_of_ints; i++) {
      tmp = 0;
      for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
        tmp = bytes[bytecnt] * sizes[i] + tmp;
        bytes[bytecnt] = tmp & 0xff;
        tmp >>= 8;
      }

      while (tmp != 0) {
        bytes[bytecnt++] = tmp & 0xff;
        tmp >>= 8;
      }
      num_of_bytes = bytecnt;
    }

    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
      num_of_bits++;
      num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}


static int sizeofints(int num_of_ints, unsigned int sizes[])
{
  alignas(64) unsigned char bytes[32];

  *(unsigned int*)bytes = sizes[0]; 
  unsigned int msb = 31 - __builtin_clz(sizes[0]); 
  unsigned int num_of_bytes = msb/8 + 1; 
  int i = 1;

#pragma GCC unroll 3 // max 4, take off our check start check
  for (; i < num_of_ints; i++) {
    unsigned int tmp = 0;
    unsigned int bytecnt = 0;

    for (; bytecnt < num_of_bytes; bytecnt++) {
      tmp += bytes[bytecnt] * sizes[i];
      bytes[bytecnt] = tmp;
      tmp >>= 8;
    }
  
    /* tmp is effectivly an overflow, write the overflow 
     * to the buffer */
    unsigned char n = 0;  
    if (tmp > 0) {
      msb = 31 - __builtin_clz(tmp); 
      n = msb/8 + 1;
    }
    
    *(unsigned int*)(bytes + bytecnt) = tmp; 
    bytecnt += n; 
    num_of_bytes = bytecnt;
  }
  
  const unsigned int num_of_bits = minbits_table[bytes[--num_of_bytes]] ; 
  const unsigned int out = num_of_bits + num_of_bytes * 8; 
  return out;
}


static int sizeofints_sse2(int num_of_ints, unsigned int sizes[]) {

  alignas(64) __m128i bytes_lane[32];
  
  /* 
   * the intial condition for each "byte chain" is a 
   * multiplicative successor e.g 
   * lane[0] = {sizes[0], sizes[0], sizes[0]*sizes[1], sizes[0]*sizes[1]*sizes[2]} 
   */
  

  return 0; 
}


int main(int argc, char *argv[])
{
  srand(time(NULL));  

  unsigned int sizes[4]; 
  for (unsigned int b=0; b<TEST_N; b++) {
    uint64_t totals = 0; 
    uint64_t times  = 0; 
    for (unsigned int i=0;i<ITER;i++) {
      sizes[0] = rand() % UINT32_MAX-1; 
      sizes[1] = rand() % UINT32_MAX-1; 
      sizes[2] = rand() % UINT32_MAX-1; 
      sizes[3] = rand() % UINT32_MAX-1; 

#if 1
      fprintf(stderr, "0x%X * 0x%X * 0x%X * 0x%X\n", 
              sizes[0], sizes[1], sizes[2], sizes[3]); 
      for (unsigned int i=0; i<4; i++) 
        print_bits_u32(sizes[i]); 
      fputc('\n', stdout); 
#endif
      
      _mm_clflush(sizes);
      uint64_t start = __rdtsc(); 
      uint64_t start_ms = now_ms(); 
      volatile unsigned int sink = sizeofints(4, sizes); 
      uint64_t end_ms   = now_ms(); 
      uint64_t end   = __rdtsc(); 

      totals += end - start; 
      times += end_ms - start_ms; 

      volatile unsigned int true = sizeofints_og(4, sizes); 
      if (sink != true) {
        fprintf(stderr, "logical error! - %u vs %u\n", sink, true);
        abort(); 
      }
    }

    fprintf(stderr, "cycles: %llu\ttime: %llu ms\n", 
            (unsigned long long)totals, (unsigned long long)times);
  }
  
  return 0; 
}
