#include<stdio.h>
#include<stdlib.h>

int main(int argc, char *argv[]){
	
	int n = atoi(argv[1]);
	int Nfourier = atoi(argv[2]);
	int FFE = atoi(argv[3]);
	int N = atoi(argv[4]);
	
	FILE *f;
	f = fopen("parameters.h", "wt");
	fprintf(f, "#ifndef COMMON_INCLUDE_FILE\n");
	fprintf(f, "#define COMMON_INCLUDE_FILE\n");
	fprintf(f, "\n");
	fprintf(f, "#define n %d //number of oscillators\n", n);
	fprintf(f, "#define Nfourier %d //fourier terms used for prc\n", Nfourier);
	fprintf(f, "#define FEE %d //first estimation of epsilon (1 -> yes, 0 -> no)\n", FFE);
	fprintf(f, "#define N %d //coarse graining of the interval 2Pi\n", N);
	fprintf(f, "\n");
	fprintf(f, "#endif\n");
	fclose(f);
	
	return 0;
	
}
