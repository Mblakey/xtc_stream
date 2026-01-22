

#include <stdio.h>



static unsigned char highest_set_bit(unsigned char x) {
  
  for (int i=7; i >= 0; i--) {
    if (x & (1 << i)) 
      return i+1;
  }
  return 0;
}


int main() {
  
  printf("unsigned char jmp_table[256] = {\n"); 
  for (unsigned short i=0; i<256; i++) {
    printf("%u,", highest_set_bit(i)); 

    if (i % 4 == 0)
      printf("\n"); 
  }
  printf("};\n"); 
  return 0; 
}
