#include<stdio.h>
#include<stdlib.h>
#include<iostream>
//#include<Eigen/Dense>
#include"../Eigen/Dense"
#include"rools.h"
#include"../parameters.h"
#define Pi 3.141592653589793238462643383

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){

	int observed = atoi(argv[1]);
	
	FILE *f;
	//read the pulses
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
			if(st[i] != 1){
				if(t == spikes[i][st[i]-2]){ //if theres a repetition
					st[i]--;
				}
			}
		}
		fclose(f);
	}
	if(st[observed] == 0){ //if theres no spikes from a unit then theres no way to tell anything
		f = fopen("EPS.txt","wt");
		for(int i = 0; i < n; ++i){
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
	//make scatter plots for codes (for preliminary eps estimation)
	int code_length = N; //just first spike counts
	double **avg = new double*[n];
	int **avg_ct = new int*[n];
	for(int i = 0; i < n; ++i){
		avg[i] = new double[code_length];
		avg_ct[i] = new int[code_length];
	}
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < code_length; ++j){
			avg[i][j] = 0;
			avg_ct[i][j] = 0;
		}
	}
	for(int i = 0; i < n; ++i){
		if(i != observed){
			for(int in = 0; in < st[observed]-1-1; ++in){ //loop over intervals
				int code = int(floor((intervals[i][in][1]/intervals[i][in][int(round(intervals[i][in][0]))+1])*N)); //just first spike counts
				if(code < code_length){ //should always be the case
					avg[i][code] += intervals[i][in][int(intervals[i][in][0])+1];
					avg_ct[i][code]++;
				}
			}
		}
	}
	//and with stdev
	double *eps_r = new double[n];
	//first calculate avg
	double *avgavg = new double[n];
	for(int i = 0; i < n; ++i){
		if(i != observed){
			avgavg[i] = 0;
			int sst = 0;
			for(int j = 0; j < code_length-1; ++j){ //you take all but the last, because the last may be half full
				if(avg_ct[i][j] != 0){
					avgavg[i] += avg[i][j]/avg_ct[i][j];
					sst++;
				}
			}
			if(sst > 0) avgavg[i] /= sst;
		}
	}
	//and the std
	for(int i = 0; i < n; ++i){
		if(i != observed){
			eps_r[i] = 0;
			int sst = 0;
			for(int j = 0; j < code_length; ++j){
				if(avg[i][j] != 0){
					eps_r[i] += pow(avg[i][j]/avg_ct[i][j]-avgavg[i],2);
					sst++;
				}
			}
			if(sst > 0) eps_r[i] /= sst;
		}
	}
	//find the largest link
	double largest = 0;
	for(int i = 0; i < n; ++i){
		if(eps_r[i] > largest) largest = eps_r[i];
	}
	
	//file print
	f = fopen("EPS.txt","wt");
	for(int i = 0; i < n; ++i){
		if(st_cap[i] == 0) fprintf(f, "%lf\n", 0.0/0.0); //if there are no pulses then its not true that there is no connections, its juts non aplicable
		else{
			if(FEE) fprintf(f, "%lf\n", sqrt(eps_r[i])/sqrt(largest));
			else fprintf(f, "%lf\n", 1.);
		}
	}
	fclose(f);
	
	return 0;
	
}

