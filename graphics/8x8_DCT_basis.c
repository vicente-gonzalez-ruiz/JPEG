#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DCT.h"

int main(int argc, char *argv[]) {
  double samples[64];   /* vector de muestras */
  unsigned char *image;

  //samples = (double *)malloc(8*8*sizeof(double));
  image = (unsigned char *)calloc(72*72,sizeof(unsigned char));

  {
    int x, y;
    for(y=0; y<8; y++) {
      for(x=0; x<8; x++) {
	fprintf(stderr,"%d %d ",x,y);
	{
	  int i;
	  for(i=0; i<64; i++) {
	    samples[i] = 0.0;
	  }
	  samples[y*8+x] = 64.0;
	  ifct2d(samples, 8, 8);
	  {
	    int i,j;
	    unsigned char *ptr = image+((y*9)*9*8+(x*9));
	    for(i=0; i<8; i++) {
	      for(j=0; j<8; j++) {
		int x = floor(samples[i*8+j])+128;
		if(x<0) x = 0;
		else if(x>255) x = 255;
		ptr[i*72+j] = x;
		fprintf(stderr,"%d ",x);
	      }
	    }
	    fprintf(stderr,"\n");
	  }
	}
      }
    }
  }
  fwrite(image, sizeof(unsigned char), 72*72, stdout);
  return 0;
}
