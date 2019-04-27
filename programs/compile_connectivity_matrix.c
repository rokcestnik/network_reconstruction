#include<stdio.h>
#include<math.h>
#include"../parameters.h"
#define Pi 3.141592653589793238462643383

double PRC(double *prc, double x){
	double res = prc[1];
	for(int p = 0; p < Nfourier; ++p){
		res += prc[1+2*p+1]*sin((p+1)*x);
		res += prc[1+2*p+2]*cos((p+1)*x);
	}
	return res;
}

double normalization(double *prc){
	double res = 0;
	for(double phi = 0; phi < 2*Pi; phi = phi+0.005){
		res += fabs(PRC(prc,phi));
	}
	return res;
}

int main(void){
	
	FILE *f;
	//just clear the files
	f = fopen("reconstruction/connectivity.txt","wt");
	fclose(f);
	f = fopen("reconstruction/PRCs.txt","wt");
	fclose(f);
	f = fopen("reconstruction/frequencies.txt","wt");
	fclose(f);
	
	//loop over nodes
	for(int i = 0; i < n; ++i){
		//read prc
		double *prc = new double[1+1+2*Nfourier];
		char buffer[1024];
		snprintf(buffer, sizeof(buffer), "reconstruction/PRC%d.txt", i);
		f = fopen(buffer,"rt");
		for(int p = 0; p < 1+1+2*Nfourier; ++p){
			fscanf(f, "%lf\n", &(prc[p]));
		}
		fclose(f);
		//calculate PRC and EPS normalization
		double c = normalization(prc);
		//read eps
		double *eps = new double[n];
		snprintf(buffer, sizeof(buffer), "reconstruction/EPS%d.txt", i);
		f = fopen(buffer,"rt");
		for(int j = 0; j < n; ++j){
			fscanf(f, "%lf\n", &(eps[j]));
		}
		fclose(f);
		//and now print
		f = fopen("reconstruction/connectivity.txt","at");
		for(int j = 0; j < n; ++j){
			fprintf(f, "%lf", eps[j]*c);
			if(j == n-1) fprintf(f, "\n");
			else fprintf(f, ", ");
		}
		fclose(f);
		f = fopen("reconstruction/PRCs.txt","at");
		for(int p = 1; p < 1+1+2*Nfourier; ++p){
			fprintf(f, "%lf", prc[p]/c);
			if(p == 1+2*Nfourier) fprintf(f, "\n");
			else fprintf(f, ", ");
		}
		fclose(f);
		f = fopen("reconstruction/frequencies.txt","at");
		fprintf(f, "%lf\n", prc[0]);
		fclose(f);
		//delete individual EPS and PRC files
		snprintf(buffer, sizeof(buffer), "reconstruction/EPS%d.txt", i);
		remove(buffer);
		snprintf(buffer, sizeof(buffer), "reconstruction/PRC%d.txt", i);
		remove(buffer);
	}
	
	
	return 0;
	
}

