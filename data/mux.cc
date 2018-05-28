#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "args.hh"

#define N_BYTES_1 1
#define N_BYTES_2 1

int main(int argc, char *argv[]) {
  args args(argc, argv);
  if(args.exist_arg((char *)"-help")) {
    cerr << "Multiplex 2 data sequences\n";
    cerr << " second input data seq. < first input data seq. > output data seq.\n";
    cerr << " [-n_bytes_1 <chunk size of the first source> (default " << N_BYTES_1 << ")\n";
    cerr << " [-n_bytes_2 <chunk size of the second source> (default " << N_BYTES_2 << ")\n";
  } else {

    FILE *input_file_1 = stdin;
    FILE *input_file_2 = fopen(argv[1],"rb");
    if(!input_file_2) {
      cerr << argv[0] << ": unable to open \"" << argv[1] << "\"\n";
      return 1;
    }
    FILE *output_file = stdout;

    int n_bytes_1 = N_BYTES_1;
    if(args.exist_arg((char *)"-n_bytes_1")) n_bytes_1 = atoi(args.get_arg());
    cerr << argv[0] << ": number of bytes of the first source = " << n_bytes_1 << '\n';

    int n_bytes_2 = N_BYTES_2;
    if(args.exist_arg((char *)"-n_bytes_2")) n_bytes_2 = atoi(args.get_arg());
    cerr << argv[0] << ": number of bytes of the second source = " << n_bytes_2 << '\n';
    
    char *buffer_1 = (char *)malloc(n_bytes_1);
    if(!buffer_1) {
      cerr << argv[0] << ": unable to allocate " << n_bytes_1 << " bytes of memory\n";
      return 1;
    }
    char *buffer_2 = (char *)malloc(n_bytes_2);
    if(!buffer_2) {
      cerr << argv[0] << ": unable to allocate " << n_bytes_2 << " bytes of memory\n";
      return 1;
    }
    for(;;) {
      fread(buffer_1, sizeof(char), n_bytes_1, input_file_1);
      if(feof(input_file_1)) break;
      fread(buffer_2, sizeof(char), n_bytes_2, input_file_2);
      if(feof(input_file_2)) break;
      fwrite(buffer_1, sizeof(char), n_bytes_1, output_file);
      fwrite(buffer_2, sizeof(char), n_bytes_2, output_file);
    }

    free(buffer_1);
    free(buffer_2);
  }
  return 0;
}
