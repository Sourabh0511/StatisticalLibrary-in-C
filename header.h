#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#define epsilon 0.00001
//#define MAX 200

double mean(double *a,int n);
void sort(double *a,int low,int high);
int partition(double *a,int low,int high);
void swap(double *a,double *b);
double median_sorted(double *a,int n);
double almost_equal(double a,double b);
double mode(double a[],int n);
double iqr(double *a,int n); 
double max(double *a,int n);
double min(double *a,int n);
double median_unsorted(double *a,int n);
//double kth_smallest(double *a,int k,int l,int r);
double quickselect(double *A, int left, int right, int k);
int partition2(double *A,int left,int right);
double summary_statistics(double *a,int n);
double copy(double *a,double *b,int n);
//int kth_select(double *a,int s,int e,int k);
double variance(double *a,int n);
double standard_deviation_sample(double *a,int n);
double test_f();
double p_norm(const double z_value);
double chisqr(int Dof, double Cv);
static double igf(double S, double Z);
int z_test_onetail(double X,double mu,double sd,double alpha,int n)	;
int z_test_twotail(double X,double mu,double sd,double alpha,int n);
double z_ci(double alpha,double mean,double stdev,int n,double *lci,double *ucl);
double prob_to_zsc(double prob);
double t_value(const double x,const double mean,const double standard_deviation,int n);
int t_test(double t_value,double alpha,int n);
int chi_sqr_goodness_of_fit(int *observ_freq,int *expected_freq,double alpha,int n);
double pearson_r(double *x,double *y,int n);
double r_squared(double *x,double *y,int n);
int least_square_fit(double *x,double *y,int n,double *m,double *b);
