/*
 *  Michael Blakey 
 *  GROMACS XTC reverse engineering
 *  for streamed decompression
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include <string.h>
#include <time.h>

#include <curl/curl.h>

#include "xdrfile.h"
#include "xdrfile_xtc.h"

size_t block_size;
bool opt_verbose; 
bool opt_timing; 
bool opt_force_stream; 
char path_type; 
const char *path; 


#define MAGIC 1995

#define PATH_TYPE_FPTR 0
#define PATH_TYPE_URL  1
#define PATH_IPV4      2
#define PATH_IPV6      3

#define BYTES_PER_XDR_UNIT 4
static char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};

static const int magicints[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
  80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
  1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
  16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031,
  131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
  832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
  4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216
};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


enum state_op {
  HEADER_READ_MAGIC,
  HEADER_READ_NATOMS,
  HEADER_READ_STEP,
  HEADER_READ_TIME,
  COORD_READ_BOX,
  COORD_READ_SIZE,
  COORD_CHECK_NATOM,
  COORD_READ_PRECISION,
  COORD_READ_MIN_INT,
  COORD_READ_MAX_INT,
  COORD_READ_SMALLIDX,
  COORD_READ_BYTE_LENGTH,
  COORD_SETUP_OPAQUE_XDRS,
  COORD_READ_OPAQUE_XDRS,
  COORD_READ_OPAQUE_CRUD,
  COORD_READ_FRAME,
} state;  


union punned_int32_t {
  uint8_t barr[4];
  int32_t _int; 
}; 


int *decomp_buffer;
size_t decomp_buffer_size; 

int *buf1; 
int *buf2; 
char *cp; 
size_t frames  = 0; 
size_t idx_ptr = 0;

int smallidx; 
int cnt; 
int lsize; 

int32_t xtc_natoms; 
float xtc_time; 
float xtc_prec; 
int32_t xtc_step; 

unsigned int rndup;
char crud[BYTES_PER_XDR_UNIT]; 

size_t lread  = 0; 
union punned_int32_t lp; 

int minint[3]; 
int maxint[3];  
float box[DIM*DIM]; 
rvec *pos = NULL; 


/* decompile the compressed xtc file directly from a data path */
static int decompile_xtcfile(struct XDRFILE *fp) 
{
  if (opt_verbose) fprintf(stderr, "=== reading xtc file...\n"); 

  clock_t start; 

  if (opt_timing) 
    start = clock(); 
  
  int natoms; 
  int step; 
  float time; 
  matrix box; 
  float prec = 1000; 
  /*
  float *lx = &box[0][0]; 
  float *ly = &box[1][1]; 
  float *lz = &box[2][2]; 
  */

  int xdr_status = xtc_header(fp, &natoms, &step, &time, 0x0);
  if (xdr_status == exdrOK){
    if (opt_verbose)
      fprintf(stderr, "=== %d atoms in xtc\n", natoms); 
  }
  else fprintf(stderr, "Error: decompile_xtcfile() -- read atoms returned non exdrOK status\n"); 

  rvec *pos = (rvec*)malloc(natoms*3*sizeof(float));
  
  /* Read frames of an open xtc file */
  xdr_seek(fp, SEEK_SET, 0); 
  while (read_xtc(fp, natoms, &step, &time, box, pos, &prec) == exdrOK) {
    frames++; 
    fprintf(stdout, "%u  %u  %f\n", natoms, step, time); 
    for (unsigned int i=0; i<natoms; i++)
      fprintf(stdout, "%8.3f %8.3f %8.3f\n",pos[i][0], pos[i][1], pos[i][2]);
  }
  
  if (opt_verbose) fprintf(stderr, "=== %lu frames read\n", frames); 
  free(pos);

  if (opt_timing) {
    clock_t end = clock(); 
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC; 

    int seconds = (int)elapsed; 
    int mseconds = (int)((elapsed - seconds) * 1000.0); 

    fprintf(stderr, "%d.%03d seconds\n", seconds, mseconds); 
  }
  return 0; 
}


static unsigned int fill_bitsize_arrays(unsigned int *bitsizeint, 
                                      unsigned int *sizeint, 
                                      int *minint, 
                                      int *maxint) 
{
  unsigned int bitsize;  
  sizeint[0] = maxint[0] - minint[0]+1;
  sizeint[1] = maxint[1] - minint[1]+1;
  sizeint[2] = maxint[2] - minint[2]+1;
  if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
    bitsizeint[0] = sizeofint(sizeint[0]);
    bitsizeint[1] = sizeofint(sizeint[1]);
    bitsizeint[2] = sizeofint(sizeint[2]);
    bitsize = 0; /* flag the use of large sizes */
  }
  else 
    bitsize = sizeofints(3,sizeint); 
  
  return bitsize; 
}


static void decompress_opaque_block()
{
  int flag; 
  int tmp;
  int is_smaller; 
  int prevcoord[3];
  unsigned int sizesmall[3]; 
  unsigned int bitsizeint[3], sizeint[3]; 

  unsigned int bitsize = fill_bitsize_arrays(bitsizeint, sizeint, minint, maxint); 

  tmp = smallidx-1;
  tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
  int smaller = magicints[tmp]/2;
  int smallnum = magicints[smallidx]/2;

  sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx]; 
  buf2[0] = buf2[1] = buf2[2] = 0;

  int *lip   = buf1;
	float *lfp = pos[0]; 
  float inv_precision = 1.0/xtc_prec;
  
  int i = 0; 
  int run = 0;
  while (i < lsize) {
	  int *thiscoord = (int*)(lip) + i*3;
    if (bitsize == 0) {
      thiscoord[0] = decodebits(buf2, bitsizeint[0]);
      thiscoord[1] = decodebits(buf2, bitsizeint[1]);
      thiscoord[2] = decodebits(buf2, bitsizeint[2]);
    }
    else 
      decodeints(buf2, 3, bitsize, sizeint, thiscoord);

    i++;
    thiscoord[0] += minint[0];
    thiscoord[1] += minint[1];
    thiscoord[2] += minint[2];

    prevcoord[0] = thiscoord[0];
    prevcoord[1] = thiscoord[1];
    prevcoord[2] = thiscoord[2];

    flag = decodebits(buf2, 1);
    is_smaller = 0;
    if (flag == 1) {
      run = decodebits(buf2, 5);
      is_smaller = run % 3;
      run -= is_smaller;
      is_smaller--;
    }
    if (run > 0) {
      thiscoord += 3;
      for (int k = 0; k < run; k+=3) {
        decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
        i++;
        thiscoord[0] += prevcoord[0] - smallnum;
        thiscoord[1] += prevcoord[1] - smallnum;
        thiscoord[2] += prevcoord[2] - smallnum;
        if (k == 0) {
          /* interchange first with second atom for better
           * compression of water molecules
           */
          tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
          prevcoord[0] = tmp;
          tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
          prevcoord[1] = tmp;
          tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
          prevcoord[2] = tmp;
          *lfp++ = prevcoord[0] * inv_precision;
          *lfp++ = prevcoord[1] * inv_precision;
          *lfp++ = prevcoord[2] * inv_precision;
        } 
        else {
          prevcoord[0] = thiscoord[0];
          prevcoord[1] = thiscoord[1];
          prevcoord[2] = thiscoord[2];
        }
        *lfp++ = thiscoord[0] * inv_precision;
        *lfp++ = thiscoord[1] * inv_precision;
        *lfp++ = thiscoord[2] * inv_precision;
      }
    }
    else {
      *lfp++ = thiscoord[0] * inv_precision;
      *lfp++ = thiscoord[1] * inv_precision;
      *lfp++ = thiscoord[2] * inv_precision;
    }

    smallidx += is_smaller;
    if (is_smaller < 0) {
      smallnum = smaller;
      if (smallidx > FIRSTIDX) 
        smaller = magicints[smallidx - 1] /2;
      else 
        smaller = 0;
    }
    else if (is_smaller > 0) {
      smaller = smallnum;
      smallnum = magicints[smallidx] / 2;
    }
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
  }

  for (unsigned int i=0; i<lsize; i++)
    fprintf(stdout, "%8.3f %8.3f %8.3f\n",pos[i][0], pos[i][1], pos[i][2]);
}


#define CHUNK_SUCCESS  0
#define CHUNK_ERROR    1
#define CHUNK_CONTINUE 2


static int decompile_xtc_chunk(uint8_t *blk, size_t gbytes)
{
  size_t gread = 0; 
  while (gread < gbytes) {

    switch (state) {

      case HEADER_READ_MAGIC:
jmp_header_read_magic:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;

        if (xdr_ntohl(lp._int) != MAGIC) {
          fprintf(stderr, "Error: decompile_blocked_xtc() -- did not read magic\n");  
          return CHUNK_ERROR; 
        }
        state = HEADER_READ_NATOMS; 
        lread = 0; 

      case HEADER_READ_NATOMS:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;

        xtc_natoms = xdr_ntohl(lp._int); 
        if (pos == NULL) 
          pos = (rvec*)malloc(xtc_natoms*3*sizeof(float));
        state = HEADER_READ_STEP; 
        lread = 0; 

      case HEADER_READ_STEP:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;

        xtc_step = xdr_ntohl(lp._int); 
        state = HEADER_READ_TIME; 
        lread = 0; 
        
      case HEADER_READ_TIME:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;
  
        // TODO: possible bug here on different architectures
        lp._int = xdr_ntohl(lp._int); 
        *(int *)&xtc_time = lp._int;

        fprintf(stdout, "%u  %u  %f\n", xtc_natoms, xtc_step, xtc_time); 

        state = COORD_READ_BOX; 
        lread   = 0; 
        idx_ptr = 0; 
      
      /* read the full box coordinates */
      case COORD_READ_BOX: 
        while (gread<gbytes && idx_ptr < DIM*DIM) {
          lp.barr[lread++] = blk[gread++];
          if (lread==4) {
            lp._int = xdr_ntohl(lp._int); 
            *(int32_t*)&box[idx_ptr++] = (int32_t)lp._int;
            lread = 0; 
          }
        }

        if (idx_ptr < DIM*DIM) 
          return CHUNK_CONTINUE;

        state = COORD_READ_SIZE; 
        lread = 0; 
        idx_ptr = 0; 

      case COORD_READ_SIZE: 
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;
        
        lsize = xdr_ntohl(lp._int);
        if (xtc_natoms < lsize) {
          fprintf(stderr, "Error: decompile_blocked_xtc() -- requested %u coords, file contains %u\n", lsize, xtc_natoms);  
          return CHUNK_ERROR; 
        }

        if (lsize*3*2 > decomp_buffer_size) {
          decomp_buffer_size = lsize*2*3; 
          decomp_buffer =  (int*)realloc(decomp_buffer, decomp_buffer_size*sizeof(int)); 
        }
        buf1 = decomp_buffer; 
        buf2 = &decomp_buffer[lsize*3];

        /* buf2[0-2] are special and do not contain actual data */
        buf2[0] = buf2[1] = buf2[2] = 0;

        state = COORD_CHECK_NATOM; 
        lread = 0; 
      
      case COORD_CHECK_NATOM:
        /* 
         * need to check this with Josh, and unit test small
         * files where this is true. e.g what are we expecting to see
         *
         *
         *  if(*size<=9)
         *    {
         *    return xdrfile_read_float(ptr,size3,xfp)/3;
         *    * return number of coords, not floats *
         *  }
         */
        if (xtc_natoms <= 9) {
          fprintf(stderr ,"[TMP ABORT] - natoms < 9, see source code for comments\n");
          abort(); 
          return CHUNK_ERROR;
        }
        state = COORD_READ_PRECISION;  
        lread = 0; 

      case COORD_READ_PRECISION:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;

        lp._int = xdr_ntohl(lp._int); 
        *(int *)&xtc_prec = lp._int;

        state = COORD_READ_MIN_INT;  
        idx_ptr = 0; 
        lread   = 0; 

      case COORD_READ_MIN_INT:
        while (gread<gbytes && idx_ptr < 3) {
          lp.barr[lread++] = blk[gread++];
          if (lread==4) {
            minint[idx_ptr++] = xdr_ntohl(lp._int); 
            lread = 0; 
          }
        }

        if (idx_ptr < 3) 
          return CHUNK_CONTINUE;

        state = COORD_READ_MAX_INT; 
        idx_ptr = 0; 
        lread   = 0; 

      case COORD_READ_MAX_INT:
        while (gread<gbytes && idx_ptr < 3) {
          lp.barr[lread++] = blk[gread++];
          if (lread==4) {
            maxint[idx_ptr++] = xdr_ntohl(lp._int); 
            lread = 0; 
          }
        }

        if (idx_ptr < 3) 
          return CHUNK_CONTINUE;

        state = COORD_READ_SMALLIDX; 
        lread = 0; 

      case COORD_READ_SMALLIDX:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
        if (lread != 4) 
          return CHUNK_CONTINUE;

        smallidx = xdr_ntohl(lp._int); 
        state = COORD_READ_BYTE_LENGTH; 
        lread = 0; 

	    /* buf2[0] holds the length in bytes */
      case COORD_READ_BYTE_LENGTH:
        while (gread<gbytes && lread<4) 
          lp.barr[lread++] = blk[gread++];
       if (lread != 4) 
          return CHUNK_CONTINUE;

        buf2[0] = xdr_ntohl(lp._int); 
        state = COORD_SETUP_OPAQUE_XDRS;
        lread = 0; 

      case COORD_SETUP_OPAQUE_XDRS:
        cnt = buf2[0]; 
	      rndup = cnt % BYTES_PER_XDR_UNIT;
        if (rndup > 0)
          rndup = BYTES_PER_XDR_UNIT - rndup;

        cp = (char*)&(buf2[3]); 
        state = COORD_READ_OPAQUE_XDRS;
        idx_ptr = 0; 
        lread   = 0; 
    
      /* 
       * streamed version of 
       * xdrfile_read_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp) 
       */
      case COORD_READ_OPAQUE_XDRS:
        while (gread<gbytes && idx_ptr < cnt) 
          cp[idx_ptr++] = blk[gread++]; 
        if (idx_ptr != cnt) 
          return CHUNK_CONTINUE;

        state = COORD_READ_OPAQUE_CRUD; 
        lread   = 0; 
        idx_ptr = 0; 

			/* return xdr_getbytes (xdrs, (char *)crud, rndup); */
      case COORD_READ_OPAQUE_CRUD:
        if (rndup == 0) {
          state = COORD_READ_FRAME;
          lread = 0; 
        }
        else {
          while (gread<gbytes && idx_ptr < rndup) 
            crud[idx_ptr++] = blk[gread++]; 

          if (idx_ptr != rndup) 
            return CHUNK_CONTINUE;

          state = COORD_READ_FRAME;
          lread = 0; 
        }

      case COORD_READ_FRAME:
        decompress_opaque_block(); 
        state = HEADER_READ_MAGIC; 
        lread   = 0; 
        idx_ptr = 0; 
        frames++; 
        goto jmp_header_read_magic;
    }
  }

  return CHUNK_SUCCESS; 
}


static int decompile_chunked_xtc_fp(FILE *fp, uint8_t *blk, size_t blk_size) 
{
  if (opt_verbose) fprintf(stderr, "=== streaming xtc file statefully...\n"); 
  
  frames = 0; 
  state = HEADER_READ_MAGIC; 
  decomp_buffer_size = 16;
  decomp_buffer = (int*)malloc(sizeof(int)*decomp_buffer_size); 
  
  size_t gbytes; 
  while ((gbytes = fread(blk, 1, blk_size, fp))!=0) {
    if (decompile_xtc_chunk(blk, gbytes) != CHUNK_CONTINUE)
      break; 
  }

  if (pos) 
    free(pos); 
  free(decomp_buffer); 
  
  if (opt_verbose) fprintf(stderr, "=== %lu frames read\n", frames); 
  return 0; 
}



static size_t libcurl_write_callback(void *ptr, size_t size, size_t nmemb, void *data)
{
  if (size*nmemb != 0) {
    if (decompile_xtc_chunk(ptr, size*nmemb) == CHUNK_ERROR)
      return 0; 
  }
  return (size_t)(size * nmemb);
}


static int decompile_chunked_xtc_curl(CURL *eh, uint8_t *blk, size_t blk_size) {
  if (opt_verbose) fprintf(stderr, "=== streaming xtc file statefully over download...\n"); 
  
  state = HEADER_READ_MAGIC; 
  decomp_buffer_size = 16;
  decomp_buffer = (int*)malloc(sizeof(int)*decomp_buffer_size); 

  curl_easy_setopt(eh, CURLOPT_WRITEFUNCTION, libcurl_write_callback); 
  curl_easy_setopt(eh, CURLOPT_URL, path);
  
  int res = curl_easy_perform(eh); 

  if (pos) 
    free(pos); 
  free(decomp_buffer); 

  if (res != CURLE_OK) return 1; 
  
  if (opt_verbose) fprintf(stderr, "=== %lu frames read\n", frames); 

  curl_off_t val;
  res = curl_easy_getinfo(eh, CURLINFO_SIZE_DOWNLOAD_T, &val);
  if ((CURLE_OK == res) && (val > 0))
    fprintf(stderr, "=== downloaded %lu bytes\n", (unsigned long)val); 

  res = curl_easy_getinfo(eh, CURLINFO_TOTAL_TIME_T, &val);
  if ((CURLE_OK == res) && (val > 0)) {
    fprintf(stderr, "=== download time: %lu.%06lu sec\n",
           (unsigned long)(val / 1000000), (unsigned long)(val % 1000000));
  }

  /* check for average download speed */
  res = curl_easy_getinfo(eh, CURLINFO_SPEED_DOWNLOAD_T, &val);
  if ((CURLE_OK == res) && (val > 0)) {
    fprintf(stderr, "=== avg download speed: %lu kbyte/sec\n",
           (unsigned long)(val / 1024));
  }
  
  return res; 
}


static void display_usage() {
  fprintf(stderr, 
      "usage:\n"
      "  xtc-decompile [options] [<path> | -u <url>]\n"
      "options:\n"
      //"  -b|--block <size>\tset the thread block size (default 16KiB)\n"
      "  -t|--timing\ttime the decompression\n"
      "  -u|--url\tstream the xtc file from a valid file url\n"
      "  -v|--verbose\tverbose logging to stderr\n"
      ); 
  exit(1); 
}


static void process_cml(int argc, char **argv) {
  int j = 0; 
  
  opt_verbose      = false; 
  opt_force_stream = false; 
  opt_timing       = false; 

  block_size = 16384; 
  path_type  = PATH_TYPE_FPTR; 
  path  = (const char*)NULL; 

  for (int i=1; i<argc; i++) {
    const char *ptr = argv[i]; 

    if (ptr[0] == '-' && ptr[1]) switch(ptr[1]) {
      case 'b':
        if (++i == argc || !(block_size = atoi(argv[i]))) {
          fprintf(stderr, "Error: process_cml -- block size option must be followed with an int\n"); 
          display_usage(); 
        }
        break; 

      case 't': opt_timing = true; break; 
      case 'u': path_type = PATH_TYPE_URL; break; 
      case 'v': opt_verbose++; break; 

      case '-':
        if (strcmp(ptr, "--force-stream") == 0)
          opt_force_stream = true;
        else if (strcmp(ptr, "--verbose") == 0)
          opt_verbose++; 
        else if (strcmp(ptr, "--timing") == 0)
          opt_timing = true; 
        else if (strcmp(ptr, "--url") == 0)
          path_type = PATH_TYPE_URL; 
        else display_usage(); 
        break; 

      default: 
        fprintf(stderr, "Error: process_cml -- unknown option - %s\n", ptr); 
        display_usage(); 
    }
    else switch (j++) {
      case 0: path = ptr; break; 
      default: display_usage(); 
    }
  }

  if (j < 1) {
    fprintf(stderr, "Error: process_cml -- not enough args\n"); 
    display_usage(); 
  }
}




int main(int argc, char **argv) {
  process_cml(argc, argv); 
 
  switch (path_type) {

    case PATH_TYPE_FPTR: 
      if (opt_force_stream) {
        FILE *ifp = fopen(path, "rb"); 
        if (ifp == NULL) {
          fprintf(stderr, "Error: main() -- could not open xtc file at %s\n", path); 
          return 1; 
        }

        uint8_t *blk = (uint8_t*)malloc(block_size); 
        decompile_chunked_xtc_fp(ifp, blk, block_size); 
        
        free(blk); 
        fclose(ifp); 
        return 0; 
      }
      else {
        XDRFILE *fp = xdrfile_open(path, "rb"); 
        if (fp == NULL) {
          fprintf(stderr, "Error: main() -- could not open xtc file at %s\n", path); 
          return 1; 
        }
        decompile_xtcfile(fp); 
        xdrfile_close(fp);
        return 0; 
      }

    case PATH_TYPE_URL: 
      curl_global_init(CURL_GLOBAL_ALL); 
      CURL *eh = curl_easy_init(); 
      uint8_t *blk = (uint8_t*)malloc(block_size); 

      decompile_chunked_xtc_curl(eh, blk, block_size); 

      free(blk); 
      curl_easy_cleanup(eh);
      curl_global_cleanup(); 
      return 0; 
  }

  return 0; 
}


