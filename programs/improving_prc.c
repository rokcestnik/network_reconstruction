#include<stdio.h>
#include<iostream>
//#include<Eigen/Dense>
#include"../Eigen/Dense"
#include"rools.h"
#include"../parameters.h"
#define Pi 3.141592653589793238462643383

using namespace std;
using namespace Eigen;

double PRC(double *prc, double x){
	double res = prc[1];
	for(int p = 0; p < Nfourier; ++p){
		res += prc[1+2*p+1]*sin((p+1)*x);
		res += prc[1+2*p+2]*cos((p+1)*x);
	}
	return res;
}

int main(int argc, char *argv[]){

	int observed = atoi(argv[1]);
	
	//first just read the prc
	double *prc = new double[1+1+2*Nfourier];
	FILE *f;
	f = fopen("PRC.txt","rt");
	for(int p = 0; p < 1+1+2*Nfourier; ++p){
		fscanf(f, "%lf\n", &(prc[p]));
		if(prc[p] != prc[p]) prc[p] = 0; //if its nan that I write 0 so it doesnt error
	}
	fclose(f);
	//eps also
	double *eps = new double[n];
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
		for(int i = 0; i < n; ++i) stimuli[i] = 0;
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
	//and now to go through the intervals, and at the same time fill A and b
	VectorXf r(1+1+2*Nfourier); //result
	if(st[observed]-1 > 2*Nfourier+1+1){ //if there are enough interspike intervals for an overdetermined system
		MatrixXf A(st[observed]-1-1, 1+1+2*Nfourier); //sepravi vzamem vse intervale, pa omege ne fitam
		VectorXf b(st[observed]-1-1);
		for(int in = 0; in < st[observed]-1-1; ++in){
			//first just count 'em (sej ves, c pa to)
			int spikes_ct = 0;
			for(int i = 0; i < n; ++i){
				if(i != observed){
					for(int local_st = 1; local_st <= intervals[i][in][0]; ++local_st){
						spikes_ct++;
					}
				}
			}
			double **sspikes = new double*[spikes_ct];
			for(int s = 0; s < spikes_ct; ++s) sspikes[s] = new double[3]; //0 -> which channel, 1 -> time, 2 -> phase
			//and now again
			spikes_ct = 0;
			for(int i = 0; i < n; ++i){
				if(i != observed){
					for(int local_st = 1; local_st <= intervals[i][in][0]; ++local_st){
						sspikes[spikes_ct][0] = i;
						sspikes[spikes_ct][1] = intervals[i][in][local_st];
						spikes_ct++;
					}
				}
			}
			merge_sort_by(sspikes,spikes_ct,3,1);
			//now lets determine the phases (from prc and eps)
			if(spikes_ct > 0){ //this only makes sense if there are any spikes in the interval
				//first rescaling for 2Pi
				double rescaling = 0;
				//first one by hand
				rescaling += sspikes[0][1]*prc[0] + eps[int(round(sspikes[0][0]))]*PRC(prc, rescaling+sspikes[0][1]*prc[0]);
				for(int s = 1; s < spikes_ct; ++s){ //and now from the second one onward
					double dphase = (sspikes[s][1]-sspikes[s-1][1])*prc[0] + eps[int(round(sspikes[s][0]))]*PRC(prc, rescaling+(sspikes[s][1]-sspikes[s-1][1])*prc[0]);
					rescaling += dphase;
				}
				if(observed != 0) rescaling += (intervals[0][in][int(round(intervals[0][in][0]))+1]-sspikes[spikes_ct-1][1])*prc[0]; //the last part between the last spike and the end of the interval I add manually
				else rescaling += (intervals[1][in][int(round(intervals[1][in][0]))+1]-sspikes[spikes_ct-1][1])*prc[0]; //the last part between the last spike and the end of the interval I add manually
				//and now we go one by one and assign phases
				double phase = 0;
				//first by hand again
				phase += sspikes[0][1]*prc[0];
				sspikes[0][2] = phase/rescaling*(2*Pi);
				phase += eps[int(round(sspikes[0][0]))]*PRC(prc, phase);
				for(int s = 1; s < spikes_ct; ++s){
					phase += (sspikes[s][1]-sspikes[s-1][1])*prc[0];
					sspikes[s][2] = phase/rescaling*(2*Pi);
					phase += eps[int(round(sspikes[s][0]))]*PRC(prc, phase);
				}
			}
			//now lets fill A
			for(int p = 0; p < 1+1+2*Nfourier; ++p){
				A(in,p) = 0;
			}
			if(observed != 0) A(in, 0) = intervals[0][in][int(round(intervals[0][in][0]))+1]; //T
			else A(in, 0) = intervals[1][in][int(round(intervals[1][in][0]))+1]; //T
			for(int s = 0; s < spikes_ct; ++s){
				double epss = eps[int(round(sspikes[s][0]))];
				A(in, 1) += epss;
				for(int p = 1; p <= Nfourier; ++p){
					A(in, 2*p) += epss*sin(p*sspikes[s][2]);
					A(in, 2*p+1) += epss*cos(p*sspikes[s][2]);
				}
			}
			//and lets fill b
			b(in) = 2*Pi;
		}
		//and minimize
		MatrixXf ATAA = A.transpose()*A;
		ATAA = ATAA.inverse();
		ATAA = ATAA*A.transpose();
		r = ATAA*b;
	}
	else{
		printf("\t\tPRC: not enough intervals for an overdetermined system\n");
		for(int i = 0; i < 1+1+2*Nfourier; ++i){
			r(i) = 0.0/0.0; //it should be clear that its not aplicable
		}
	}
	
	//file print
	f = fopen("PRC.txt","wt");
	for(int i = 0; i < 1+1+2*Nfourier; ++i){
		fprintf(f, "%lf\n", r(i));
	}
	fclose(f);
	
	return 0;
	
}

