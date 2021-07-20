#include "header.h"

int main()
{
	int n;
	scanf("%d",&n);
	double *a=(double*)malloc(sizeof(double)*n);
	int i;
	time_t t;
	//printf("%lf\n",prob_to_zsc(0.0777));
	/*double lci;
	double uci;
	z_ci(0.80,50,2,100,&lci,&uci);
	printf("%lf %lf\n",lci,uci);*/
	//printf("%d\n",t_test(2.75,0.05,10));
	// double x[10]={1, 2, 4,  5,  10, 20};
	// double y[10]={4, 6, 12, 15, 34, 68};
	// double m=0.0;
	// double b=0.0;
	// least_square_fit(x,y,6,&m,&b);
	// printf("%lf\n",m);
	// printf("%lf\n",b);
	// printf("%lf\n",pearson_r(x,y,6));
	double b[11]={1.0,2.0,2.0,2.0,3.0,4.0,4.0,4.0,1.0,2.0};
	printf("%lf\n",mode(b,11));
	srand((unsigned) time(&t));
	 for(i=0;i<n;i++)
	 {
			
	 	int r=rand()%1000000;
	 	a[i]=(double)(r/10000.0);
	 }
	 for(i=0;i<n;i++)
	 {
	 	printf("%lf ",a[i]);
	 }
	
	 printf("\n");
	 int mid=n/2;
	//printf("%lf\n",quickselect(a,0,n-1,mid+1));
	//printf("%lf\n",median_unsorted(a,n));
	sort(a,0,n-1);
	summary_statistics(a,n);
	printf("p value : %lf \n",p_norm(-2.1834));
	//printf("%lf",kth_smallest(a,3,0,n-1));
	/*printf("%lf\n",variance(a,n));
	printf("%lf",standard_deviation_sample(a,n));
	printf("\n");*/

	/*
	printf("p value for chi2 %lf\n",chisqr(20,12.0));
	printf("Z test\n");*/

	//printf("%d\n",z_test_onetail(4.5,5.4,2.7,0.05,80));
	//printf("%d\n",z_test_twotail(4.5,5.4,2.7,0.05,80));
	//printf("test : %lf",test_f());
	//sort(a,0,n-1);
	//

	for(i=0;i<n;i++)
	{
		printf("%lf ",a[i]);
	}
	printf("\n");
	//printf("%lf ",max(a,n));
	//printf("%lf ",min(a,n));
	//printf("Mode = %lf",mode(a,n));
	return 0;
}
