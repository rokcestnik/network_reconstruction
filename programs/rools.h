/*

wait(time_t delay)
	waits delay seconds

clear_file(const char *fileName)
	clears the contents of the file
	
copy_file(const char *fileName1, const char *fileName2)
	copies file1 to file2
	
write_on_file_text(const char *fileName, const char *contents)
	appends *contents to the file
	
write_on_file_int(const char *fileName, int v)
	appends v to the file

write_on_file_double(const char *fileName, double v)
	appends *contents to the file
	
file_length(const char *fileName)
	returns the number of characters in a file

file_lines(const char *fileName)
	returns the number of lines in a file
		ex:
			printf("number of lines = %d\n", file_lines("dat.txt"));

max_index(double *tab, int n)
	returns the index of the largest element in table *tab with n elements

min_index(double *tab, int n)
	returns the index of the smallest element in table *tab with n elements
	
max_element(double *tab, int n)
	returns the largest element of the table *tab with n elements
	
min_element(double *tab, int n)
	returns the smallest element of the table *tab with n elements
	
average(double *tab, int n)
	returns the average of the table *tab with n elements
	
in_list(double *list, int n, double element)
	returns 1 if element is in the list, and 0 if it is not. *list needs to be SORTED!!

histogram(double *samples, int n, int *hist, int NObins)
	makes a histogram of *samples (size n) with NObins bins

linear_fit(double **points, int n, double *A, double *B, double *error)
	fits a line Ax+b through n points
		ex:
			double A, B, er;
			linear_fit(p, n, &A, &B, &er);

bubble_sort(double *tab, int n)
	bubble sorts the table *tab with n elements
	
merge_sort(double *tab, int n)
	merge sorts the table *tab with n elements
	
merge_sort_memory(double **tab, int n)
	the same as merge_sort(), but uses temporary hard memory
		ex:
			int n = 100;
			double *tab;
			tab = new double[n];
			for(int i = 0; i < n; ++i) tab[i] = rrand();
			merge_sort_memory(&tab, n);
	
merge_sort_by(double **tab, int n, int m, int element)
	merge sorts the table **tab with respect to "element" (sorting by either the first, second, third,... element in the row). 
	**tab has n rows and m elements in each row
	
merge_sort_by_memory(double ***tab, int n, int m, int element)
	the same as merge_sort_by(), but uses temporary hard memory (see merge_sort_memory() for analogy)
	
autocorrelation(double *sample, int n, double taumin, double taumax, double factor, const char *fileName)
	calculates the autocorrelation function of "sample" (size n) for tau ranging from "taumin" to "taumax" by "factor" and writes it to a file

cross_correlation(double *sample1, double *sample2, int n, double taumin, double taumax, double factor, const char *fileName)
	calculates the cross-correlation function of "sample1" and "sample2" (both size n) for tau ranging from "taumin" to "taumax" by "factor" and writes it to a file

linear_interpolation(double x)
	linearly interpolates the function between given data points (data point have to be ordered in an ascending trend)
		ex:
			linear_interpolation_init("data.txt"); //the data has to be in the form "%lf, %lf\n"
			printf("%lf\n", linear_interpolation(0.618));

linear_interpolation_3D(double x, double y)
	linearly interpolates the function between given data points (data point have to be ordered in an ascending trend in both x and y)
		ex:
			linear_interpolation_3D_init("data.txt"); //the data has to be in the form "%lf, %lf, %lf\n"
			printf("%lf\n", linear_interpolation_3D(0.618, 1.618));
			
*/

#ifndef ROOLS_H
#define ROOLS_H

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

void beep(){
	printf("%c", 7);
}

void print_int(int x){
	printf("%d\n", x);
}

void print_double(double x){
	printf("%lf\n", x);
}

void wait(time_t delay){
	time_t time0, time1;
	time(&time0);
	do time(&time1);
	while(time1-time0 < delay);
}

void clear_file(const char *fileName){
	FILE *f = fopen(fileName, "wt");
	fclose(f);
}

void copy_file(const char *fileName1, const char *fileName2){
	FILE *f1 = fopen(fileName1, "rt");
	FILE *f2 = fopen(fileName2, "wt");
	char c;
	while(true){
		if(fscanf(f1, "%c", &c) == -1) break;
		fprintf(f2, "%c", c);
	}
	fclose(f1);
	fclose(f2);
}

void write_on_file_text(const char *fileName, const char *contents){
	FILE *f = fopen(fileName, "at");
	fputs(contents, f);
	fclose(f);
}

void write_on_file_int(const char *fileName, int v){
	FILE *f = fopen(fileName, "at");
	fprintf(f, "%d", v);
	fclose(f);
}

void write_on_file_double(const char *fileName, double v){
	FILE *f = fopen(fileName, "at");
	fprintf(f, "%lf", v);
	fclose(f);
}

int file_length(const char *fileName){
	FILE *f;
	f = fopen(fileName, "rt");
	int st = 0;
	char tmp;
	while(fscanf(f, "%c", &tmp) != EOF)	st++;
	fclose(f);
	return st;
}

int file_lines(const char *fileName){
	FILE *f;
	f = fopen(fileName, "rt");
	int st = 0;
	char tmp;
	while(fscanf(f, "%c", &tmp) != EOF) if(tmp == '\n') st++;
	fclose(f);
	return st;
}

int max_index(double *tab, int n){
	int index = 0;
	double m = tab[0];
	for(int i = 1; i < n; ++i){
		if(tab[i] > m){
			m = tab[i];
			index = i;
		}
	}
	return index;
}

int min_index(double *tab, int n){
	int index = 0;
	double m = tab[0];
	for(int i = 1; i < n; ++i){
		if(tab[i] < m){
			m = tab[i];
			index = i;
		}
	}
	return index;
}

double max_element(double *tab, int n){
	double m = tab[0];
	for(int i = 1; i < n; ++i){
		if(tab[i] > m) m = tab[i];
	}
	return m;
}

double min_element(double *tab, int n){
	double m = tab[0];
	for(int i = 1; i < n; ++i){
		if(tab[i] < m) m = tab[i];
	}
	return m;
}

double average(double *tab, int n){
	double sum = 0;
	for(int i = 0; i < n; ++i){
		sum += tab[i];
	}
	return sum/n;
}

int in_list(double *list, int n, double element){
	//bisection
	double dindex = n/2.0;
	double range = n/4.0;
	do{
		int index = int(floor(dindex));
		if(list[index] == element) return 1;
		if(list[index] < element) dindex += range;
		if(list[index] > element) dindex -= range;
		range /= 2;
	}while(range > 0.25);
	return 0;
}

void histogram(double *samples, int n, int *hist, int NObins){
	//determining the bin
	double max = -1000000000;
	double min = 1000000000;
	for(int i = 0; i < n; ++i){
		if(samples[i] > max) max = samples[i];
		if(samples[i] < min) min = samples[i];
	}
	double bin = (1.01*max-min)/NObins;
	//filling the histogram
	for(int i = 0; i < NObins; ++i) hist[i] = 0;
	for(int i = 0; i < n; ++i){
		hist[int((samples[i]-min)/bin)]++;
	}
}

void linear_fit(double **points, int nn, double *AA, double *BB, double *error){
	double x = 0;
	double y = 0;
	double xy = 0;
	double x2 = 0;
	for(int i = 0; i < nn; ++i){
		x += points[i][0];
		y += points[i][1];
		xy += points[i][0]*points[i][1];
		x2 += pow(points[i][0],2);
	}
	*AA = (nn*xy/x-y)/(nn*x2/x-x);
	*BB = (xy-*AA*x2)/x;
	*error = 0;
	for(int i = 0; i < nn; ++i) *error += pow(points[i][1]-(*AA*points[i][0]+*BB),2);
}

void bubble_sort(double *tab, int n){
	for(int i = 0; i < n-1; ++i){
		for(int j = 0; j < n-i-1; ++j){
			if(tab[j] > tab[j+1]){
				double tmp = tab[j];
				tab[j] = tab[j+1];
				tab[j+1] = tmp;
			}
		}
	}
}

void merge_sort(double *tab, int n){
	if(n > 1){
		double *tab1 = new double[n/2];
		double *tab2 = new double[n-(n/2)];
		for(int i = 0; i < n/2; ++i) tab1[i] = tab[i];
		for(int i = 0; i < n-(n/2); ++i) tab2[i] = tab[(n/2)+i];
		//rekurzija
		merge_sort(tab1, n/2);
		merge_sort(tab2, n-(n/2));
		//dokler ENA OD OBEH ne pride do konca
		int i = 0, j = 0;
		while(i < n/2 && j < n-(n/2)){
			if(tab1[i] < tab2[j]){
				tab[i+j] = tab1[i];
				++i;
			}
			else{
				tab[i+j] = tab2[j];
				++j;
			}
		}
		//potem pesebej spraznim drugo
		if(i < n/2){
			for(int k = i; k < n/2; ++k)
			tab[k+j] = tab1[k];
		}
		if(j < n-(n/2)){
			for(int k = j; k < n-(n/2); ++k)
			tab[i+k] = tab2[k];
		}
		delete[] tab1;
		delete[] tab2;
	}
}

void merge_sort_memory(double **tab, int n){
	if(n > 1){
		double *tab1, *tab2;
		if(n > 100){ //ce je velik zapisem na file, da pol loh zbrisem iz rama
			FILE *f;
			f = fopen("tmp_memory.txt","wt");
			for(int i = 0; i < n; ++i) fprintf(f,"%.12lf ", (*tab)[i]);
			fclose(f);
			delete[] (*tab);
			tab1 = new double[n/2];
			tab2 = new double[n-(n/2)];
			f = fopen("tmp_memory.txt","rt");
			for(int i = 0; i < n/2; ++i) fscanf(f,"%lf ", &tab1[i]);
			for(int i = 0; i < n-(n/2); ++i) fscanf(f,"%lf ", &tab2[i]);
			fclose(f);
			remove("tmp_memory.txt");
		}
		else{ //ce je premejhn kr v ramu
			tab1 = new double[n/2];
			tab2 = new double[n-(n/2)];
			for(int i = 0; i < n/2; ++i) tab1[i] = (*tab)[i];
			for(int i = 0; i < n-(n/2); ++i) tab2[i] = (*tab)[(n/2)+i];
		}
		//rekurzija
		merge_sort_memory(&tab1, n/2);
		merge_sort_memory(&tab2, n-(n/2));
		//zdj pa realokejtam (samo ce je n > 100)
		if(n > 100){
			(*tab) = new double[n];
		}
		//dokler ENA OD OBEH ne pride do konca
		int i = 0, j = 0;
		while(i < n/2 && j < n-(n/2)){
			if(tab1[i] < tab2[j]){
				(*tab)[i+j] = tab1[i];
				++i;
			}
			else{
				(*tab)[i+j] = tab2[j];
				++j;
			}
		}
		//potem pesebej spraznim drugo
		if(i < n/2){
			for(int k = i; k < n/2; ++k)
			(*tab)[k+j] = tab1[k];
		}
		if(j < n-(n/2)){
			for(int k = j; k < n-(n/2); ++k)
			(*tab)[i+k] = tab2[k];
		}
		delete[] tab1;
		delete[] tab2;
	}
}

void merge_sort_by(double **tab, int n, int m, int element){
	if(n > 1){
		double **tab1 = new double *[n/2];
		double **tab2 = new double *[n-(n/2)];
		for(int i = 0; i < n/2; ++i) tab1[i] = new double[m];
		for(int i = 0; i < n-(n/2); ++i) tab2[i] = new double[m];
		for(int i = 0; i < n/2; ++i) for(int j = 0; j < m; ++j) tab1[i][j] = tab[i][j];
		for(int i = 0; i < n-(n/2); ++i) for(int j = 0; j < m; ++j) tab2[i][j] = tab[(n/2)+i][j];
		//rekurzija
		merge_sort_by(tab1, n/2, m, element);
		merge_sort_by(tab2, n-(n/2), m, element);
		//dokler ENA OD OBEH ne pride do konca
		int i = 0, j = 0;
		while(i < n/2 && j < n-(n/2)){
			if(tab1[i][element] < tab2[j][element]){
				for(int mm = 0; mm < m; ++mm) tab[i+j][mm] = tab1[i][mm];
				++i;
			}
			else{
				for(int mm = 0; mm < m; ++mm) tab[i+j][mm] = tab2[j][mm];
				++j;
			}
		}
		//potem pesebej spraznim drugo
		if(i < n/2){
			for(int k = i; k < n/2; ++k){
				for(int mm = 0; mm < m; ++mm) tab[k+j][mm] = tab1[k][mm];
			}
		}
		if(j < n-(n/2)){
			for(int k = j; k < n-(n/2); ++k){
				for(int mm = 0; mm < m; ++mm) tab[i+k][mm] = tab2[k][mm];
			}
		}
		for(int i = 0; i < n/2; ++i) delete[] tab1[i];
		for(int i = 0; i < n-(n/2); ++i) delete[] tab2[i];
		delete[] tab1;
		delete[] tab2;
	}
}

void merge_sort_by_memory(double ***tab, int n, int m, int element){
	if(n > 1){
		double **tab1, **tab2;
		if(n > 100){ //ce je velik zapisem na file, da pol loh zbrisem iz rama
			FILE *f;
			f = fopen("tmp_memory.txt","wt");
			for(int i = 0; i < n; ++i){
				for(int j = 0; j < m; ++j){
					fprintf(f,".12%lf ", (*tab)[i][j]);
				}
			}
			fclose(f);
			for(int i = 0; i < n; ++i) delete[] (*tab)[i];
			delete[] (*tab);
			tab1 = new double*[n/2];
			for(int i = 0; i < n/2; ++i) tab1[i] = new double[m];
			tab2 = new double*[n-(n/2)];
			for(int i = 0; i < n-(n/2); ++i) tab2[i] = new double[m];
			f = fopen("tmp_memory.txt","rt");
			for(int i = 0; i < n/2; ++i){
				for(int j = 0; j < m; ++j){
					fscanf(f,"%lf ", &tab1[i][j]);
				}
			}
			for(int i = 0; i < n-(n/2); ++i){
				for(int j = 0; j < m; ++j){
					fscanf(f,"%lf ", &tab2[i][j]);
				}
			}
			fclose(f);
			remove("tmp_memory.txt");
		}
		else{ //ce je premejhn kr v ramu
			tab1 = new double*[n/2];
			for(int i = 0; i < n/2; ++i) tab1[i] = new double[m];
			tab2 = new double*[n-(n/2)];
			for(int i = 0; i < n-(n/2); ++i) tab2[i] = new double[m];
			for(int i = 0; i < n/2; ++i){
				for(int j = 0; j < m; ++j){
					tab1[i][j] = (*tab)[i][j];
				}
			}
			for(int i = 0; i < n-(n/2); ++i){
				for(int j = 0; j < m; ++j){
					tab2[i][j] = (*tab)[(n/2)+i][j];
				}
			}
		}
		//rekurzija
		merge_sort_by_memory(&tab1, n/2, m, element);
		merge_sort_by_memory(&tab2, n-(n/2), m, element);
		//zdj pa realokejtam (samo ce je n > 100)
		if(n > 100){
			(*tab) = new double*[n];
			for(int i = 0; i < n; ++i) (*tab)[i] = new double[m];
		}
		//dokler ENA OD OBEH ne pride do konca
		int i = 0, j = 0;
		while(i < n/2 && j < n-(n/2)){
			if(tab1[i][element] < tab2[j][element]){
				for(int mm = 0; mm < m; ++mm) (*tab)[i+j][mm] = tab1[i][mm];
				++i;
			}
			else{
				for(int mm = 0; mm < m; ++mm) (*tab)[i+j][mm] = tab2[j][mm];
				++j;
			}
		}
		//potem pesebej spraznim drugo
		if(i < n/2){
			for(int k = i; k < n/2; ++k){
				for(int mm = 0; mm < m; ++mm) (*tab)[k+j][mm] = tab1[k][mm];
			}
		}
		if(j < n-(n/2)){
			for(int k = j; k < n-(n/2); ++k){
				for(int mm = 0; mm < m; ++mm) (*tab)[i+k][mm] = tab2[k][mm];
			}
		}
		for(int i = 0; i < n/2; ++i) delete[] tab1[i];
		for(int i = 0; i < n-(n/2); ++i) delete[] tab2[i];
		delete[] tab1;
		delete[] tab2;
	}
}

void autocorrelation(double *sample, int n, double taumin, double taumax, double factor, const char *fileName){
	double avg = 0; //average
	for(int i = 0; i < n; ++i) avg += sample[i];
	avg /= n;
	double var = 0; //variance
	for(int i = 0; i < n; ++i) var += pow(sample[i]-avg,2);
	FILE *f;
	//the actual correlation
	for(double dtau = taumin; dtau < taumax; dtau *= factor){
		int tau = int(round(dtau));
		double correlation = 0;
		for(int i = 0; i < n; ++i) correlation += (sample[i]-avg)*(sample[(i+tau)%n]-avg);
		correlation /= var;
		//file print
		f = fopen(fileName,"at");
		fprintf(f,"%d,%lf\n", tau, correlation);
		fclose(f);
	}
}

void cross_correlation(double *sample1, double *sample2, int n, double taumin, double taumax, double factor, const char *fileName){
	double avg1 = 0; //average1
	for(int i = 0; i < n; ++i) avg1 += sample1[i];
	avg1 /= n;
	double avg2 = 0; //average2
	for(int i = 0; i < n; ++i) avg2 += sample2[i];
	avg2 /= n;
	double var1 = 0; //variance1
	for(int i = 0; i < n; ++i) var1 += pow(sample1[i]-avg1,2);
	double var2 = 0; //variance2
	for(int i = 0; i < n; ++i) var2 += pow(sample2[i]-avg2,2);
	FILE *f;
	//the actual correlation
	for(double dtau = taumin; dtau < taumax; dtau *= factor){
		int tau = int(round(dtau));
		double correlation = 0;
		for(int i = 0; i < n; ++i) correlation += (sample1[i]-avg1)*(sample2[(i+tau)%n]-avg2);
		correlation /= pow(var1*var2,0.5);
		//file print
		f = fopen(fileName,"at");
		fprintf(f,"%d,%lf\n", tau, correlation);
		fclose(f);
	}
}

double linear_interpolation_table[100000][2];
int linear_interpolation_table_length;
void linear_interpolation_init(const char *filename){
	FILE *f;
	f = fopen(filename, "rt");
	int st = 0;
	while(fscanf(f, "%lf, %lf\n", &linear_interpolation_table[st][0], &linear_interpolation_table[st][1]) == 2) st++;
	linear_interpolation_table_length = st;
	fclose(f);
}
double linear_interpolation(double x){
	if(x < linear_interpolation_table[0][0] || x > linear_interpolation_table[linear_interpolation_table_length-1][0]) return 0;
	//bisection
	double index = linear_interpolation_table_length/2.;
	double interval = linear_interpolation_table_length/4.;
	while(!(linear_interpolation_table[int(index)][0] <= x && x <= linear_interpolation_table[int(index)+1][0])){
		if(x < linear_interpolation_table[int(index)][0]) index -= interval;
		else index += interval;
		interval /= 2;
	}
	double koeficient = (linear_interpolation_table[int(index)+1][1]-linear_interpolation_table[int(index)][1])/(linear_interpolation_table[int(index)+1][0]-linear_interpolation_table[int(index)][0]);
	return linear_interpolation_table[int(index)][1]+koeficient*(x-linear_interpolation_table[int(index)][0]);
}

double linear_interpolation_3D_table[200][200][3];
int linear_interpolation_3D_table_length_x;
int linear_interpolation_3D_table_length_y;
void linear_interpolation_3D_init(const char *filename){
	FILE *f;
	f = fopen(filename, "rt");
	//we scan the whole file into a table
	double tab[40000][3];
	int st = 0;
	while(fscanf(f, "%lf, %lf, %lf\n", &tab[st][0], &tab[st][1], &tab[st][2]) == 3) st++;
	fclose(f);
	//now we parse the data to the "linear_interpolation_3D_table"
	int st_x = 0;
	int st_y = 0;
	linear_interpolation_3D_table[0][0][0] = tab[0][0];
	linear_interpolation_3D_table[0][0][1] = tab[0][1];
	linear_interpolation_3D_table[0][0][2] = tab[0][2];
	//if the data is in the form:
	//											x0, y0, z
	//											x0, y1, z
	//											x0, y2, z
	//											x1, y0, z
	//											x1, y1, z
	//											x1, y2, z
	if(tab[0][0] == tab[1][0]){
		for(int i = 1; i < st; ++i){
			if(tab[i][0] == tab[i-1][0]){
				st_y++;
			}
			else{
				st_x++;
				st_y = 0;
			}
			linear_interpolation_3D_table[st_x][st_y][0] = tab[i][0];
			linear_interpolation_3D_table[st_x][st_y][1] = tab[i][1];
			linear_interpolation_3D_table[st_x][st_y][2] = tab[i][2];
		}
	}
	//if the data is in the form:
	//											x0, y0, z
	//											x1, y0, z
	//											x2, y0, z
	//											x0, y1, z
	//											x1, y1, z
	//											x2, y1, z
	if(tab[0][1] == tab[1][1]){
		for(int i = 1; i < st; ++i){
			if(tab[i][1] == tab[i-1][1]){
				st_x++;
			}
			else{
				st_y++;
				st_x = 0;
			}
			linear_interpolation_3D_table[st_x][st_y][0] = tab[i][0];
			linear_interpolation_3D_table[st_x][st_y][1] = tab[i][1];
			linear_interpolation_3D_table[st_x][st_y][2] = tab[i][2];
		}
	}
	linear_interpolation_3D_table_length_x = st_x+1;
	linear_interpolation_3D_table_length_y = st_y+1;
	//print
	/*printf("x\\y\n");
	for(int x = 0; x < linear_interpolation_3D_table_length_x; ++x){
		for(int y = 0; y < linear_interpolation_3D_table_length_y; ++y){
			printf("%d\\%d\t%lf, %lf, %lf\n", x, y, linear_interpolation_3D_table[x][y][0],linear_interpolation_3D_table[x][y][1],linear_interpolation_3D_table[x][y][2]);
		}
	}*/
}
double linear_interpolation_3D(double x, double y){
	if(x < linear_interpolation_3D_table[0][0][0] || x > linear_interpolation_3D_table[linear_interpolation_3D_table_length_x-1][0][0] || y < linear_interpolation_3D_table[0][0][1] || y > linear_interpolation_3D_table[0][linear_interpolation_3D_table_length_y-1][1]) return 0;
	//bisection x
	double index_x = linear_interpolation_3D_table_length_x/2.;
	double interval_x = linear_interpolation_3D_table_length_x/4.;
	while(!(linear_interpolation_3D_table[int(index_x)][0][0] <= x && x <= linear_interpolation_3D_table[int(index_x)+1][0][0])){
		if(x < linear_interpolation_3D_table[int(index_x)][0][0]) index_x -= interval_x;
		else index_x += interval_x;
		interval_x /= 2;
	}
	//bisection y
	double index_y = linear_interpolation_3D_table_length_y/2.;
	double interval_y = linear_interpolation_3D_table_length_y/4.;
	while(!(linear_interpolation_3D_table[0][int(index_y)][1] <= y && y <= linear_interpolation_3D_table[0][int(index_y)+1][1])){
		if(y < linear_interpolation_3D_table[0][int(index_y)][1]) index_y -= interval_y;
		else index_y += interval_y;
		interval_y /= 2;
	}
	double koeficient1 = (linear_interpolation_3D_table[int(index_x)+1][int(index_y)][2]-linear_interpolation_3D_table[int(index_x)][int(index_y)][2])/(linear_interpolation_3D_table[int(index_x)+1][int(index_y)][0]-linear_interpolation_3D_table[int(index_x)][int(index_y)][0]);
	double koeficient2 = (linear_interpolation_3D_table[int(index_x)+1][int(index_y)+1][2]-linear_interpolation_3D_table[int(index_x)][int(index_y)+1][2])/(linear_interpolation_3D_table[int(index_x)+1][int(index_y)+1][0]-linear_interpolation_3D_table[int(index_x)][int(index_y)+1][0]);
	double point1 = linear_interpolation_3D_table[int(index_x)][int(index_y)][2]+koeficient1*(x-linear_interpolation_3D_table[int(index_x)][int(index_y)][0]);
	double point2 = linear_interpolation_3D_table[int(index_x)][int(index_y)+1][2]+koeficient2*(x-linear_interpolation_3D_table[int(index_x)][int(index_y)+1][0]);
	double the_koeficient = (point2-point1)/(linear_interpolation_3D_table[int(index_x)][int(index_y)+1][1]-linear_interpolation_3D_table[int(index_x)][int(index_y)][1]);
	return point1+the_koeficient*(y-linear_interpolation_3D_table[int(index_x)][int(index_y)][1]);	
}

#endif
