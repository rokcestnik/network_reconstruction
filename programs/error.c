#include<stdio.h>
#include"rools.h"
#include"../parameters.h"
#define Pi 3.141592653589793238462643383

double PRC(double *prc, double x){
	double res = prc[0]; //in other scripts it might be prc[0], thats because I used to save frequencies at prc[0]
	for(int p = 0; p < Nfourier; ++p){
		res += prc[2*p+1]*sin((p+1)*x);
		res += prc[2*p+2]*cos((p+1)*x);
	}
	return res;
}

int main(void){
	
	FILE *f;
	//just clear the files from any previous runs
	f = fopen("reconstruction/err.txt","wt");
	fclose(f);
	f = fopen("reconstruction/rec_score.txt","wt");
	fclose(f);
	//global standard deviations of psi
	int psi_st = 0; //psi counter
	int psi0_st = 0; //psi0 counter
	double psi_std_global = 0;
	double psi0_std_global = 0;
	//first read connectivity
	double **eps = new double*[n];
	for(int i = 0; i < n; ++i) eps[i] = new double[n];
	f = fopen("reconstruction/connectivity.txt","rt");
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			fscanf(f, "%lf", &(eps[i][j]));
			if(eps[i][j] != eps[i][j]) eps[i][j] = 0; //if its nan we just put 0 so it doesnt crash
			if(j == n-1) fscanf(f, "\n");
			else fscanf(f, ", ");
		}
	}
	fclose(f);
	//read prcs
	double **prc = new double*[n];
	for(int i = 0; i < n; ++i) prc[i] = new double[1+2*Nfourier];
	f = fopen("reconstruction/PRCs.txt","rt");
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 1+2*Nfourier; ++j){
			fscanf(f, "%lf", &(prc[i][j]));
			if(prc[i][j] != prc[i][j]) prc[i][j] = 0; //if its nan just put 0 so it doesnt crash
			if(j == 2*Nfourier) fscanf(f, "\n");
			else fscanf(f, ", ");
		}
	}
	fclose(f);
	//read frequencies
	double *om = new double[n];
	f = fopen("reconstruction/frequencies.txt","rt");
	for(int i = 0; i < n; ++i){
		fscanf(f, "%lf\n", &(om[i]));
		if(om[i] != om[i]) om[i] = 0; //if its nan just put 0 so it doesnt crash
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
	//here is going to be the main loop over nodes
	for(int observed = 0; observed < n; ++observed){
		double psi_std;
		double psi0_std;
		if(st[observed] != 0){
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
						while(spikes[i][st[i]] < spikes[observed][in] && st[i] < st_cap[i]) st[i]++; //dangerous line
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
			//and now to go through the intervals and calculate phase at the end of intervals
			psi_std = 0;
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
				double psi = 0; //phase at end of interval
				if(spikes_ct > 0){ //when there are some spikes
					//first one by hand
					psi += sspikes[0][1]*om[observed] + eps[observed][int(round(sspikes[0][0]))]*PRC(prc[observed], psi+sspikes[0][1]*om[observed]);
					for(int s = 1; s < spikes_ct; ++s){ //and now from the second one onward
						double dphase = (sspikes[s][1]-sspikes[s-1][1])*om[observed] + eps[observed][int(round(sspikes[s][0]))]*PRC(prc[observed], psi+(sspikes[s][1]-sspikes[s-1][1])*om[observed]);
						psi += dphase;
					}
					if(observed != 0) psi += (intervals[0][in][int(round(intervals[0][in][0]))+1]-sspikes[spikes_ct-1][1])*om[observed]; //the last part between the last spike and the end of the interval I add manually
					else psi += (intervals[1][in][int(round(intervals[1][in][0]))+1]-sspikes[spikes_ct-1][1])*om[observed]; //the last part between the last spike and the end of the interval I add manually
				}
				else if(spikes_ct == 0){ //if there arent any spikes
					if(observed != 0) psi = om[observed]*intervals[0][in][1]; //mors pogledat v intervals od enga k ni observed na prvo (ki je zadno) mesto (ker 0 je kolk jih je)
					else psi = om[observed]*intervals[1][in][1];
				}
				psi_std += pow(psi-2*Pi,2);
				psi_std_global += pow(psi-2*Pi,2);
				psi_st++;
			}
			psi_std = sqrt(psi_std/(st[observed]-1-1));
			//and now also psi0
			double avg_om = 0; //average frequency
			for(int in = 0; in < st[observed]-1-1; ++in){
				if(observed != 0) avg_om += 2*Pi/intervals[0][in][int(round(intervals[0][in][0]))+1]; // 2Pi/T
				else avg_om += 2*Pi/intervals[1][in][int(round(intervals[1][in][0]))+1];
			}
			avg_om /= (st[observed]-1-1);
			//now get psi0
			double psi0;
			psi0_std = 0;
			for(int in = 0; in < st[observed]-1-1; ++in){
				if(observed != 0) psi0 = avg_om*intervals[0][in][int(round(intervals[0][in][0]))+1];
				else psi0 = avg_om*intervals[1][in][int(round(intervals[1][in][0]))+1];
				psi0_std += pow(psi0-2*Pi,2);
				psi0_std_global += pow(psi0-2*Pi,2);
				psi0_st++;
			}
			psi0_std = sqrt(psi0_std/(st[observed]-1-1));
		}
		else{
			psi_std = 0.0/0.0;
			psi0_std = 0.0/0.0;
		}
		//file print
		f = fopen("reconstruction/err.txt","at");
		fprintf(f, "%lf, %lf\n", psi0_std, psi_std);
		fclose(f);
		f = fopen("reconstruction/rec_score.txt","at");
		fprintf(f, "%lf\n", log(psi0_std/psi_std));
		fclose(f);
	}
	psi_std_global = sqrt(psi_std_global/psi_st);
	psi0_std_global = sqrt(psi0_std_global/psi0_st);
	//and for the global error
	f = fopen("reconstruction/err.txt","at");
	fprintf(f, "%lf, %lf\n", psi0_std_global, psi_std_global);
	fclose(f);
	f = fopen("reconstruction/rec_score.txt","at");
	fprintf(f, "%lf\n", log(psi0_std_global/psi_std_global));
	fclose(f);
	
	return 0;
	
}

