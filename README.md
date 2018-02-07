## phase retrieval benchmarks
### RRR: Simple Phase Retrieval
This program is deliberately minimalist so as not to obscure the structure of the algorithm.
It needs the [FFTW3 library](http://www.fftw.org).

The program is set up for solving the benchmarks described in:
> "Benchmark problems for phase retrieval", V. Elser, T.-Y. Lan & T. Bendory

To compile:
```
gcc -O2 RRR.c -lm -lfftw3 -o RRR
```

To run:
```
./RRR [datafile] [supp] [powgoal] [beta] [iterlimit] [trials] [resultfile] &
```

- datafile:	one of the benchmark datafiles (data/data100E, data/data140E, ...)
- supp:		support size = 8*N, N = 100, 140, ... is the number of atoms
- powgoal:	fractional power in support
- beta:		RRR parameter
- iterlimit:	RRR iteration limit (long int)
- trials:		number of random starts
- resultfile:	ASCII file of the iteration counts for each trial

Example:
```
./RRR data/data100E 800 .95 .5 1000 5 results100E &
```

Solutions are written to a file named sol (M x M table of floats).

Send comments, bug reports, to: ve10@cornell.edu
