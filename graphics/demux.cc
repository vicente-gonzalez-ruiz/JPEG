#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "args.hh"

#define CHUNK_SIZE 3
#define CHUNK_OFFSET 0
#define N_BYTES 1

int main(int argc, char *argv[]) {
  args args(argc, argv);
  if(args.exist_arg((char *)"-help")) {
    cerr << "Demultiplex 2 data sequences\n";
    cerr << " < multi-component input data seq. > single-component output data seq.\n";
    cerr << " [-chunk_size <multi-component chunk size in bytes>] (default " << CHUNK_SIZE << ")\n";
    cerr << " [-chunk_offset <internal chunk offset in of bytes>] (default " << CHUNK_OFFSET << ")\n";
    cerr << " [-n_bytes <number of bytes to extract per chunk>] (default " << N_BYTES << ")\n";
  } else {
    
    int chunk_size = CHUNK_SIZE;
    if(args.exist_arg((char *)"-chunk_size")) chunk_size = atoi(args.get_arg());
    cerr << argv[0] << ": chunk size = " << chunk_size << '\n';

    int chunk_offset = CHUNK_OFFSET;
    if(args.exist_arg((char *)"-chunk_offset")) chunk_offset = atoi(args.get_arg());
    cerr << argv[0] << ": internal chunk offset = " << chunk_offset << '\n';

    int n_bytes = N_BYTES;
    if(args.exist_arg((char *)"-n_bytes")) n_bytes = atoi(args.get_arg());
    cerr << argv[0] << ": number of extracted bytes/chunk = " << n_bytes << '\n';

    FILE *input_file = stdin;
    FILE *output_file = stdout;

    char *buffer = (char *)malloc(chunk_size);
    if(!buffer) {
      cerr << argv[0] << ": unable to allocate " << chunk_size << " bytes of memory\n";
      return 1;
    }
    for(;;) {
      fread(buffer, sizeof(char), chunk_size, input_file);
      if(feof(input_file)) break;
      fwrite(buffer+chunk_offset, sizeof(char), n_bytes, output_file);
    }
  }
}
