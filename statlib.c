#include "header.h"
double z_table[][10]={
	{0.5000,
	0.5040,
	0.5080,
	0.5120,
	0.5160,
	0.5199,
	0.5239,
	0.5279,
	0.5319,
	0.5359},

	{0.5398,
	0.5438,
	0.5478,
	0.5517,
	0.5557,
	0.5596,
	0.5636,
	0.5675,
	0.5714,
	0.5753},

	{0.5793,
	0.5832,
	0.5871,
	0.5910,
	0.5948,
	0.5987,
	0.6026,
	0.6064,
	0.6103,
	0.6141},

	{0.6179,
	0.6217,
	0.6255,
	0.6293,
	0.6331,
	0.6368,
	0.6406,
	0.6443,
	0.6480,
	0.6517},

	{0.6554,
	0.6591,
	0.6628,
	0.6664,
	0.6700,
	0.6736,
	0.6772,
	0.6808,
	0.6844,
	0.6879},

	{0.6915,
	0.6950,
	0.6985,
	0.7019,
	0.7054,
	0.7088,
	0.7123,
	0.7157,
	0.7190,
	0.7224},

	{0.7257,
	0.7291,
	0.7324,
	0.7357,
	0.7389,
	0.7422,
	0.7454,
	0.7486,
	0.7517,
	0.7549},

	{0.7580,
	0.7611,
	0.7642,
	0.7673,
	0.7704,
	0.7734,
	0.7764,
	0.7794,
	0.7823,
	0.7852},

	{0.7881,
	0.7910,
	0.7939,
	0.7967,
	0.7995,
	0.8023,
	0.8051,
	0.8078,
	0.8106,
	0.8133},

	{0.8159,
	0.8186,
	0.8212,
	0.8238,
	0.8264,
	0.8289,
	0.8315,
	0.8340,
	0.8365,
	0.8389},

	{0.8413,
	0.8438,
	0.8461,
	0.8485,
	0.8508,
	0.8531,
	0.8554,
	0.8577,
	0.8599,
	0.8621},

	{0.8643,
	0.8665,
	0.8686,
	0.8708,
	0.8729,
	0.8749,
	0.8770,
	0.8790,
	0.8810,
	0.8830},

	{0.8849,
	0.8869,
	0.8888,
	0.8907,
	0.8925,
	0.8944,
	0.8962,
	0.8980,
	0.8997,
	0.9015},

	{0.9032,
	0.9049,
	0.9066,
	0.9082,
	0.9099,
	0.9115,
	0.9131,
	0.9147,
	0.9162,
	0.9177},

	{0.9192,
	0.9207,
	0.9222,
	0.9236,
	0.9251,
	0.9265,
	0.9279,
	0.9292,
	0.9306,
	0.9319},

	{0.9332,
	0.9345,
	0.9357,
	0.9370,
	0.9382,
	0.9394,
	0.9406,
	0.9418,
	0.9429,
	0.9441},

	{0.9452,
	0.9463,
	0.9474,
	0.9484,
	0.9495,
	0.9505,
	0.9515,
	0.9525,
	0.9535,
	0.9545},

	{0.9554,
	0.9564,
	0.9573,
	0.9582,
	0.9591,
	0.9599,
	0.9608,
	0.9616,
	0.9625,
	0.9633},

	{0.9641,
	0.9649,
	0.9656,
	0.9664,
	0.9671,
	0.9678,
	0.9686,
	0.9693,
	0.9699,
	0.9706},

	{0.9713,
	0.9719,
	0.9726,
	0.9732,
	0.9738,
	0.9744,
	0.9750,
	0.9756,
	0.9761,
	0.9767},

	{0.9772,
	0.9778,
	0.9783,
	0.9788,
	0.9793,
	0.9798,
	0.9803,
	0.9808,
	0.9812,
	0.9817},

	{0.9821,
	0.9826,
	0.9830,
	0.9834,
	0.9838,
	0.9842,
	0.9846,
	0.9850,
	0.9854,
	0.9857},

	{0.9861,
	0.9864,
	0.9868,
	0.9871,
	0.9875,
	0.9878,
	0.9881,
	0.9884,
	0.9887,
	0.9890},

	{0.9893,
	0.9896,
	0.9898,
	0.9901,
	0.9904,
	0.9906,
	0.9909,
	0.9911,
	0.9913,
	0.9916},

	{0.9918,
	0.9920,
	0.9922,
	0.9925,
	0.9927,
	0.9929,
	0.9931,
	0.9932,
	0.9934,
	0.9936},

	{0.9938,
	0.9940,
	0.9941,
	0.9943,
	0.9945,
	0.9946,
	0.9948,
	0.9949,
	0.9951,
	0.9952},

	{0.9953,
	0.9955,
	0.9956,
	0.9957,
	0.9959,
	0.9960,
	0.9961,
	0.9962,
	0.9963,
	0.9964},

	{0.9965,
	0.9966,
	0.9967,
	0.9968,
	0.9969,
	0.9970,
	0.9971,
	0.9972,
	0.9973,
	0.9974},

	{0.9974,
	0.9975,
	0.9976,
	0.9977,
	0.9977,
	0.9978,
	0.9979,
	0.9979,
	0.9980,
	0.9981},

	{0.9981,
	0.9982,
	0.9982,
	0.9983,
	0.9984,
	0.9984,
	0.9985,
	0.9985,
	0.9986,
	0.9986},

	{0.9987,
	0.9987,
	0.9987,
	0.9988,
	0.9988,
	0.9989,
	0.9989,
	0.9989,
	0.9990,
	0.9990}

};

double t_table[][10]={{0.325, 1.000,3.078,6.314, 12.706, 31.821, 63.657, 318.309, 636.619},
{0.289,0.816,1.886,2.920,4.303,6.965,9.925,22.327,31.599},
{0.277,0.765,1.638,2.353,3.182,4.541,5.841,10.215,12.924},
{0.271,0.741,1.533,2.132,2.776,3.747,4.604,7.173,8.610},
{0.267,0.727,1.476,2.015,2.571,3.365,4.032,5.893,6.869},//5
{0.265,0.718,1.440,1.943,2.447,3.143,3.707,5.208,5.959},
{0.263,0.711,1.415,1.895,2.365,2.998,3.499,4.785,5.408},
{0.262,0.706,1.397,1.860,2.306,2.896,3.355,4.501,5.041},
{0.261,0.703,1.383,1.833,2.262,2.821,3.250,4.297,4.781},
{0.260,0.700,1.372,1.812,2.228,2.764,3.169,4.144,4.587},
{0.260,0.697,1.363,1.796,2.201,2.718,3.106,4.025,4.437},
{0.259,0.695,1.356,1.782,2.179,2.681,3.055,3.930,4.318},
{0.259,0.694,1.350,1.771,2.160,2.650,3.012,3.852,4.221},
{0.258,0.692,1.345,1.761,2.145,2.624,2.977,3.787,4.140},
{0.258,0.691,1.341,1.753,2.131,2.602,2.947,3.733,4.073},
{0.258,0.690,1.337,1.746,2.120,2.583,2.921,3.686,4.015},
{0.257,0.689,1.333,1.740,2.110,2.567,2.898,3.646,3.965},
{0.257,0.688,1.330,1.734,2.101,2.552,2.878,3.610,3.922},
{0.257,0.688,1.328,1.729,2.093,2.539,2.861,3.579,3.883},
{0.257,0.687,1.325,1.725,2.086,2.528,2.845,3.552,3.850},
{0.257,0.686,1.323,1.721,2.080,2.518,2.831,3.527,3.819},
{0.256,0.686,1.321,1.717,2.074,2.508,2.819,3.505,3.792},
{0.256,0.685,1.319,1.714,2.069,2.500,2.807,3.485,3.768},
{0.256,0.685,1.318,1.711,2.064,2.492,2.797,3.467,3.745},
{0.256,0.684,1.316,1.708,2.060,2.485,2.787,3.450,3.725},
{0.256,0.684,1.315,1.706,2.056,2.479,2.779,3.435,3.707},
{0.256,0.684,1.314,1.703,2.052,2.473,2.771,3.421,3.690},
{0.256,0.683,1.313,1.701,2.048,2.467,2.763,3.408,3.674},
{0.256,0.683,1.311,1.699,2.045,2.462,2.756,3.396,3.659},
{0.256,0.683,1.310,1.697,2.042,2.457,2.750,3.385,3.646}};

double t_table_header[]={0.40,0.25,0.10,0.05,0.025,0.01,0.005,0.001,0.0005};

double copy(double *a,double *b,int n)
{
	int i;
	for(i=0;i<n;i++){
		b[i]=a[i];
	}
}

double mean(double *a,int n){
	int i;
	double sum=0;
	for(i=0;i<n;i++){
		sum=a[i] + sum;
	}
	return sum/n;
}

void sort(double *a,int low,int high){
	if(low<high)
	{
		int p=partition(a,low,high);
		sort(a,low,p);
		sort(a,p+1,high);
	}
}

void swap(double *a,double *b){
	double temp=*a;
	*a=*b;
	*b=temp;
}

int partition(double *a,int low,int high){
	double pivot=a[low];
	int i=low -1;
	int j=high +1;
	while(1)
	{
		do{
			++i;
		}while(a[i]<pivot&&i<high);

		do{
			--j;
		}while(a[j]>pivot&&j>low);

		if(i>=j){
			return j;
		}
		swap(&a[i],&a[j]);
	}

}

double median_sorted(double *a,int n){
	int lhs=(n-1)/2;
	int rhs=n/2;
	if(n==0){
		return 0.0;
	}
	else if(lhs==rhs){
		return a[lhs];
	}
	else{
		return (a[lhs] + a[rhs])/2.0;
	}
}




double almost_equal(double a,double b)
{
    return fabs(a-b)<epsilon;

}
double mode(double a[],int n)
{
	int maxCount = 0, i, j;
	double maxValue = 0.0;

   for (i = 0; i < n; ++i) {
      int count = 0;
      
      for (j = 0; j < n; ++j) 
      {
         if (almost_equal(a[j],a[i]))
         	++count;
      }
      
      if (count > maxCount)
       {
         maxCount = count;
         maxValue = a[i];
      }
   }
   if(maxCount==1){
   	return -1;
   }
   return maxValue;
}


double iqr(double *a,int n)
{
	sort(a,0,n-1);
	double q1;
	double q3;
	int lhs=(n-1)/2;
	int rhs=n/2;
	if(lhs==rhs)
	{
		q1=a[n/4];
		q3=a[(3*n)/4];
	}
	else
	{
		q1=(a[n/4 -1]+a[n/4])/2;
		q3=(a[(3*n)/4 -1] + a[(3*n)/4])/2;
	}
	return q3-q1;
}

double max(double *a,int n)
{
	int i;
	double max=a[0];
	for(i=1;i<n;i++){
		if(a[i]>max){
			max=a[i];
		}
	}
	return max;
}

double min(double *a,int n)
{
	int i;
	double min=a[0];
	for(i=1;i<n;i++){
		if(a[i]<min){
			min=a[i];
		}
	}
	return min;
}

double median_unsorted(double *a,int n)
{
	int mid=n/2;
	if(n%2){
		return quickselect(a,0,n-1,mid+1);
	}
	else{
		return (quickselect(a,0,n-1,mid) + quickselect(a,0,n-1,mid+1))/2;
	}
}



double quickselect(double *A, int left, int right, int k){

    //p is position of pivot in the partitioned array
    int p = partition2(A, left, right);

    //k equals pivot 
    if (p == k-1){
        return A[p];
    }
    //k less than pivot
    else if (k - 1 < p){
        return quickselect(A, left, p - 1, k);
    }
    //k greater than pivot
    else{
        return quickselect(A, p + 1, right, k);
    }
}


int partition2(double *A, int left, int right){
    double pivot = A[right];
    int i = left;
    int x;
    for (x = left; x < right; x++){
        if (A[x] < pivot||almost_equal(A[x],pivot)){
            swap(&A[i], &A[x]);
            i++;
        }
    }

    swap(&A[i], &A[right]);
    return i;
}


double variance(double *a,int n)
{
	if(n==0){
		return -1;
	}
	if(n==1){
		return a[0];
	}
	double x_bar=mean(a,n);
	int i;
	double sum=0;
	for(i=0;i<n;i++)
	{
		sum=sum + (a[i]-x_bar)*(a[i]-x_bar);
	}
	double variance=sum/(n-1);
	return variance;
}

double standard_deviation_sample(double *a,int n)
{
	double var=variance(a,n);
	return sqrt(var);
}

double z_value(const double x,const double mean,const double standard_deviation)
{
	return (x-mean)/standard_deviation;
}

double summary_statistics(double *a,int n)
{
	double *b=(double*)malloc(sizeof(double)*(n+1));
	copy(a,b,n);
	sort(b,0,n-1);
	double sample_mean=mean(a,n);
	double median=median_sorted(a,n);
	double sample_mode=mode(a,n);
	double var=variance(a,n);
	double stdev=standard_deviation_sample(a,n);
	double sample_iqr=iqr(a,n);
	
	printf("------------------------------------------------\n");
	printf("Mean : %lf\n",sample_mean);
	printf("Median : %lf\n",median);
	printf("Variance : %lf\n",var);
	printf("Standard Deviation : %lf\n",stdev);
	printf("Inter Quartile Range : %lf\n",sample_iqr);
	printf("Maximum : %lf\n",b[n-1]);
	printf("Minimum : %lf\n",b[0]);
	printf("------------------------------------------------\n");

}

double p_norm(const double z_value)
{
	double z_score;
	if(z_value>0.00000){
		int z_int=round(z_value*100);
		int column=z_int%10;
		int row=z_int/10;
		if(row<=30){
			z_score=z_table[row][column];
		}
		else{
			z_score=-1.0;
		}
		
	}
	else{
		int z_int=z_value*(-100);
		int column=z_int%10;
		int row=z_int/10;
		if(row<=30){
			z_score=1.0 -z_table[row][column];
		}
		else{
			z_score=-1.0;
		}

	}
	return z_score;
}

double test_f(){
	return z_table[6][1];
}

double chisqr(int Dof, double Cv)
{
    //printf("Dof:  %i\n", Dof);
    //printf("Cv:  %f\n", Cv);
    if(Cv < 0 || Dof < 1)
    {
        return 0.0;
    }
	double K = ((double)Dof) * 0.5;		// F(x;k)= incompgamma(k/2,x/2)/tgamma(k/2)
	double X = Cv * 0.5;
	if(Dof == 2)
	{
		return exp(-1.0 * X);
	}

	double PValue = igf(K, X);
    if(isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
    {
        return 1e-14;
    }
	
	PValue /= tgamma(K);  // <math.h> implementation of the gamma function 

	return (1.0 - PValue);

}

static double igf(double S, double Z)	// incomplete gamma function.
{
	if(Z < 0.0)
	{
		return 0.0;
	}
	long double Sc = (1.0 / S);
	Sc *= powl(Z, S);
	Sc *= expl(-Z);

	long double Sum = 1.0;
	long double Nom = 1.0;
	long double Denom = 1.0;
	int I;
	for(I = 0; I < 200; I++) // 200 is the number of iterations which gives the desired precision..this value was decided through experimentation(already done by someone)
	{
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	}

	return Sum * Sc;
}

int z_test_onetail(double X,double mu,double sd,double alpha,int n)		// returns 1 if null hypothesis is rejected.Alternate is true.. mu<mu0
{
	double z_score=(X-mu)/(sd/(sqrt(n)));
	printf("%lf\n",z_score);
	double pvalue=p_norm(z_score);
	if(z_score>0.00000){
		pvalue=1.0 - pvalue;
	}
	printf("The p value:testing %lf\n",pvalue);
	if(pvalue<alpha||almost_equal(pvalue,alpha)){
		return 1;
	}
	return 0;
}

int z_test_twotail(double X,double mu,double sd,double alpha,int n)
{
	double z_score=(X-mu)/(sd/(sqrt(n)));
	double pvalue=p_norm(z_score);
	if(z_score>0.0000){
		pvalue=2.0*(1.0-pvalue);
	}
	
	//printf("The p value:testing %lf\n",pvalue);
	if(pvalue<alpha||almost_equal(pvalue,alpha)){
		return 1;
	}
	return 0;
}


double z_ci(double alpha,double mean,double stdev,int n,double *lcl,double *ucl){
	alpha=alpha + (1.0-alpha)/2;
	double z_score=prob_to_zsc(alpha);
	*lcl=mean - z_score*(stdev/sqrt(n));
	*ucl=mean + z_score*(stdev/sqrt(n));


}

double prob_to_zsc(double prob)
{
	int i=0;
	int j=0;
	int negative=0;
	double diff=1.0;
	if(prob<0.50){
		prob=1.0-prob;
		negative=1;
	}

	while(i<30 && (z_table[i][0]<prob || almost_equal(z_table[i][0],prob))){
		i++;
		
	}
	i--;
	printf("\n");
	while(j<10 && diff>fabs(z_table[i][j]-prob)){
		
		diff=fabs(z_table[i][j]-prob);
		
		j++;
	}
	if(negative){
		return -i/10.0 -(j-1)/100.0;
	}
	return i/10.0 + (j-1)/100.0;
}


double t_value(const double x,const double mean,const double standard_deviation,int n)
{
	return (x-mean)/(standard_deviation/sqrt(n));
}

int t_test(double t_value,double alpha,int n)
{
	t_value=fabs(t_value);
	int dof=n-1;
	int index=dof-1;	// because index starts from 0
	int i=0;
	while(i<9 && t_table[index][i]<t_value){
		i++;
		
	}
	if(i!=0)
	{
			double upper_bound=t_table_header[i-1];
			double lower_bound=t_table_header[i];
			printf("%lf\n",t_table_header[i-1]);
			if(almost_equal(alpha,t_table_header[i-1])||t_table_header[i-1]<alpha){
				return 1;
			}
			else{
				return 0;
			}
			//printf("%lf < P(%lf) < %lf",lower_bound,t_value,upper_bound);
	}
	else{
			
			double upper_bound=t_table_header[i];
			if(almost_equal(alpha,t_table_header[i])||t_table_header[i]<alpha){
				return 1;
			}
			else{
				return 0;
			}
			//printf("%lf < P(%lf) < %lf",lower_bound,t_value,upper_bound);
	}
	/*for(i=0;i<9;i++){
		if(t_table[index][i]>t_value||almost_equal(t_table[index][i],t_value)){
			if(i!=0){
				double lower_bound=t_table_header[i-1];
				double upper_bound=t_table_header[i];
				printf("%lf < P(%lf) < %lf",lower_bound,t_value,upper_bound);
			}
			else{
				double lower_bound=0.0;
				double upper_bound=t_table_header[i];
				printf("%lf < P(%lf) < %lf",lower_bound,t_value,upper_bound);
			}
		}
	}*/
}

int chi_sqr_goodness_of_fit(int *observ_freq,int *expected_freq,double alpha,int n)
{
	
	double chi_sqr=0.0;
	int i;
	int dof=n-1;
	for(i=0;i<n;i++){
		chi_sqr=chi_sqr + ((observ_freq[i]-expected_freq[i])*(observ_freq[i]-expected_freq[i]))/((double)expected_freq[i]);
	}
	
	double pvalue=chisqr(dof,chi_sqr);
	
	if(pvalue<alpha){
		return 1;
	}
	return 0;

}

double pearson_r(double *x,double *y,int n)	// correlation co-efficient
{
	double x_mean=mean(x,n);
	double y_mean=mean(y,n);
	double xt=0.0;
	double yt=0.0;
	double syy=0.0;
	double sxx=0.0;
	double sxy=0.0;
	int j;
	for(j=0;j<n;j++){
		xt=x[j]-x_mean;
		yt=y[j]-y_mean;
		sxx=sxx + xt*xt;
		syy=syy + yt*yt;
		sxy=sxy + xt*yt;
	}
	return sxy/(sqrt(sxx*syy));
}

double r_squared(double *x,double *y,int n)
{
	double r=pearson_r(x,y,n);
	return r*r;
}

int least_square_fit(double *x,double *y,int n,double *m,double *b) 
{
	double sx=0.0;	// sum of x
	double sx2=0.0;	//sum of x**2
	double sxy=0.0;	//sum of x*y
	double sy=0.0;	//sum of y
	double sy2=0.0;	//sum of y**2
	int i;
	for(i=0;i<n;i++){
		sx=sx + x[i];
		sx2=sx2	+ x[i]*x[i];
		sy=sy + y[i];
		sy2=sy2	+ y[i]*y[i];
		sxy=sxy + x[i]*y[i];
	}
	double denom=(n*sx2 - (sx*sx));
	if(denom==0){
		*m=0;
		*b=0;
		return 1;
	}
	*m=(n*sxy -sx*sy)/denom;
	*b=(sy*sx2 -sx*sxy)/denom;
	return 0;

}

double predict_y(double x,double m,double b)
{
	return x*m +b;
}