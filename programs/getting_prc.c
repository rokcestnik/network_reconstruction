#include<stdio.h>
#include<iostream>
//#include<Eigen/Dense>
#include"../Eigen/Dense"
#include"../parameters.h"
#define Pi 3.141592653589793238462643383

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){

	int observed = atoi(argv[1]);
	
	//first just read the eps
	double *eps = new double[n];
	FILE *f;
	f = fopen("EPS.txt","rt");
	for(int i = 0; i < n; ++i){
		fscanf(f, "%lf\n", &(eps[i]));
		if(eps[i] != eps[i]) eps[i] = 0; //if its nan that I write 0 so it doesnt error
	}
	fclose(f);
	//and now the pulses
	int *st = new int[n];
	for(int i = 0; i < n; ++i) st[i] = 0;
	for(int i = 0; i < n; ++i){
		char buffer[1024];
		snprintf(buffer, sizeof(buffer), "data/%d.txt", i);
		f = fopen(buffer,"rt");
		double t;
		while(fscanf(f, "%lf\n", &t) != EOF){
			st[i]++;
		}
		fclose(f);
	}
	//organizing the spikes
	double **spikes = new double*[n];
	for(int i = 0; i < n; ++i){
		spikes[i] = new double[st[i]];
	}
	for(int i = 0; i < n; ++i) st[i] = 0;
	for(int i = 0; i < n; ++i){
		char buffer[1024];
		snprintf(buffer, sizeof(buffer), "data/%d.txt", i);
		f = fopen(buffer,"rt");
		double t;
		while(fscanf(f, "%lf\n", &t) != EOF){
			spikes[i][st[i]] = t;
			st[i]++;
			if(st[i] != 0){
				if(t == spikes[i][st[i]-2]){ //if theres a repetition
					st[i]--;
				}
			}
		}
		fclose(f);
	}
	if(st[observed] == 0){ //if theres no spikes from a unit then theres no way to tell anything
		f = fopen("PRC.txt","wt");
		for(int i = 0; i < 1+1+2*Nfourier; ++i){
		fprintf(f, "%lf\n", 0.0/0.0);
		}
		fclose(f);
		return 0;
	}
	//and now make intervals
	double ***intervals = new double**[n];
	for(int i = 0; i < n; ++i) if(i != observed) intervals[i] = new double*[st[observed]-1]; //as many as there are intervals from the observed one
	int *st_cap = new int[n]; //these counters will keep the total number of spikes and will not be reset
	for(int i = 0; i < n; ++i) st_cap[i] = st[i];
	for(int i = 0; i < n; ++i) if(i != observed) st[i] = 0; //the counter of the observed 'st[observed]' is not reset!
	for(int in = 0; in < st[observed]-1-1; ++in){ //loop over intervals
		int *stimuli = new int[n];
		for(int i = 0; i < n; ++i) if(i != observed) stimuli[i] = 0;
		for(int i = 0; i < n; ++i){
			if(i != observed){
				while(spikes[i][st[i]] < spikes[observed][in] && st[i] < st_cap[i]) st[i]++; 
				while(spikes[i][st[i]] < spikes[observed][in+1] && st[i] < st_cap[i]){
					stimuli[i]++;
					st[i]++;
				}
			}
		}
		//allocate so much space
		for(int i = 0; i < n; ++i){
			if(i != observed) intervals[i][in] = new double[stimuli[i]+1+1]; //[0] is the number of spikes, [1] = t1, [2] = t2, [3] = t3... [-1] = T
		}
		//and fill it in
		for(int i = 0; i < n; ++i){
			if(i != observed){
				intervals[i][in][0] = stimuli[i];
				intervals[i][in][stimuli[i]+1] = spikes[observed][in+1]-spikes[observed][in]; //T
				for(int j = 1; j < stimuli[i]+1; ++j){
					intervals[i][in][j] = spikes[i][st[i]-stimuli[i]+(j-1)]-spikes[observed][in]; //t1,t2,t3...
				}
			}
		}
	}
	//lets make the matrix
	VectorXf r(1+1+2*Nfourier); //result
	if(st[observed]-1 > 2*Nfourier+1+1){ //if there are enough interspike intervals for an overdetermined system
		MatrixXf A(st[observed]-1, 1+1+2*Nfourier);
		for(int in = 0; in < st[observed]-1-1; ++in){ //loop over intervals
			if(observed != 0) A(in,0) = intervals[0][in][int(intervals[0][in][0])+1]; //T
			else A(in,0) = intervals[1][in][int(intervals[1][in][0])+1]; //T
			A(in,1) = 0;
			for(int i = 0; i < n; ++i){
				if(i != observed) A(in,1) += intervals[i][in][0]*eps[i]; //stimuli weighted with the link strength
			}
			for(int f = 0; f < Nfourier; ++f){
				double sums = 0;
				double sumc = 0;
				for(int i = 0; i < n; ++i){
					if(i != observed){
						for(int stim = 1; stim < intervals[i][in][0]+1; ++stim){
							double phase = (intervals[i][in][stim]/intervals[i][in][int(intervals[i][in][0])+1])*2*Pi; //t/T*2*Pi
							sums += sin((f+1)*phase)*eps[i];
							sumc += cos((f+1)*phase)*eps[i];
						}
					}
				}
				A(in,2+2*f) = sums;
				A(in,2+2*f+1) = sumc;
			}
		}
		//and now the result
		VectorXf b(st[observed]-1);
		for(int i = 0; i < st[observed]-1; ++i) b(i) = 2*Pi;
		MatrixXf ATAA = A.transpose()*A;
		ATAA = ATAA.inverse();
		ATAA = ATAA*A.transpose();
		r = ATAA*b;
	}
	else{
		for(int i = 0; i < 1+1+2*Nfourier; ++i){
			r(i) = 0.0/0.0; //it should be clear that its not aplicable
		}
	}
	
	//file print
	f = fopen("PRC.txt","wt");
	for(int i = 0; i < 1+1+2*Nfourier; ++i){
		fprintf(f,"%lf\n", r(i));
	}
	fclose(f);
	
	return 0;
	
}

