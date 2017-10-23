#include "func.h"
#include "util.h"

void func0(double *weights, double *arrayX, double *arrayY, int xr, int yr, int n)
{
  //pull out the division 
        int i;
	double weight = 1/((double)(n));
#pragma omp parallel for
	for(i = 0; i < n; i++){
	  weights[i] = weight;
	
		arrayX[i] = xr;
		arrayY[i] = yr;
	}
}

void func1(int *seed, int *array, double *arrayX, double *arrayY,
			double *probability, double *objxy, int *index,
			int Ones, int iter, int X, int Y, int Z, int n)
{
  //  omp_set_num_threads(12);
	int i, j;
   	int index_X, index_Y;
	int max_size = X*Y*Z;
	int roundX, roundY;
	double ones = (double)Ones;
	double acc;
	int temp;
#pragma omp parallel for  
   	for(i = 0; i < n; i++){
   		arrayX[i] += 1 + 5*rand2(seed, i);
   		arrayY[i] += -2 + 2*rand2(seed, i);
   	}
#pragma omp parallel for private(index_X, index_Y,j,acc,roundX,roundY,temp)
   	for(i = 0; i<n; i++){
	  //#pragma omp parallel for
	  roundX = round(arrayX[i]);
	  roundY = round(arrayY[i]);
	  temp = i*Ones;
	  for(j = 0; j < Ones; j++){
   			index_X = roundX + objxy[j*2 + 1];
   			index_Y = roundY + objxy[j*2];
   			index[temp + j] = fabs(index_X*Y*Z + index_Y*Z + iter);
   			if(index[temp + j] >= max_size)
   				index[temp + j] = 0;
   		}
	  //delete probability[i] = 0;
		acc = 0;
		//#pragma omp parallel for reduction (+:acc)
   		for(j = 0; j < Ones; j++) {
		  //  probability[i] += (pow((array[index[i*Ones + j]] - 100),2) -
		  //pow((array[index[i*Ones + j]]-228),2))/50.0;
		  acc += (pow((array[index[temp + j]] - 100),2) - 
		  	  pow((array[index[temp + j]]-228),2))/50.0;
		}
		acc = acc /ones; 
   		//probability[i] = probability[i]/((double) Ones);
		probability[i] = acc;
	}
}

void func2(double *weights, double *probability, int n)
{
	int i;
	double sumWeights=0;
	double temp;
#pragma omp parallel for private(temp) reduction(+:sumWeights)
	for(i = 0; i < n; i++){
	  temp = weights[i] * exp(probability[i]);
	  weights[i] = temp;
	  sumWeights += temp;

	}

#pragma omp parallel for
	for(i = 0; i < n; i++)
   		weights[i] = weights[i]/sumWeights;
}

void func3(double *arrayX, double *arrayY, double *weights, double *x_e, double *y_e, int n)
{
	double estimate_x=0.0;
	double estimate_y=0.0;
    int i;

#pragma omp parallel for reduction(+:estimate_x, estimate_y)
	for(i = 0; i < n; i++){
   		estimate_x += arrayX[i] * weights[i];
   		estimate_y += arrayY[i] * weights[i];
   	}

	*x_e = estimate_x;
	*y_e = estimate_y;

}

void func4(double *u, double u1, int n)
{
	int i;
	double N = (double)n;
#pragma omp parallel for
	for(i = 0; i < n; i++){
  		u[i] = u1 + i/N;
   	}
}

void func5(double *x_j, double *y_j, double *arrayX, double *arrayY, double *weights, double *cfd, double *u, int n)
{
	int i, j;
	double N = (double)n;
#pragma omp parallel for private(i)
	for(j = 0; j < n; j++){
   		//i = findIndex(cfd, n, u[j]);
   		i = findIndexBin(cfd, 0, n, u[j]);
   		if(i == -1)
   			i = n-1;
   		x_j[j] = arrayX[i];
   		y_j[j] = arrayY[i];

   	}
#pragma omp parallel for
	for(i = 0; i < n; i++){
		arrayX[i] = x_j[i];
		arrayY[i] = y_j[i];
		weights[i] = 1/N;
	}
}
