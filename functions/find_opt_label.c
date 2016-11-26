#include "mex.h" /* Always include this */
#define X prhs[0]


void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
	int nrhs, const mxArray *prhs[]) /* Input variables */
{
	double *Psort, *pOut, *Sum, *Num;
//	mxArray *Sum;
	int M;
	int i,j;
	double r1, r2, Den, Sum1, Sum2;


	//if (nrhs != 2 || mxDOUBLE_CLASS != mxGetClassID(X)) {
	//	mexErrMsgTxt("Number of inputs should be two!");
	//}

	Psort = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]); /* Get the dimensions of Input */
	Sum = mxGetPr(prhs[1]);
	Num = mxGetPr(prhs[2]);


	plhs[0] = mxCreateDoubleMatrix(M-1, 1, mxREAL);
	pOut = mxGetPr(plhs[0]); /* Get the pointer to the data of B */

	Sum1 = 0;
	for (i = 0; i != M-1; i++)
	{
		r1 = (((double)i) + 1) /  (double)M;
		r2 = (double)(M - i - 1) /(double)M;
//		mexPrintf("%d\n,",M);

		Sum1 = Sum1 + Psort[i];
		Sum2 = Sum[0] - Sum1;
		Den = Sum1 / (double)(i + 1) - Sum2 / (double)(M - i - 1);
		Den = Den*Den;

		if (r1 < 0.5)
		{
			pOut[i] = Num[0] / (Den*r1) - r2;
		}
		else
		{
			pOut[i] = Num[0] / (Den*r2) - r1;
		}
	}

	//for i = 1:N - 1
	//	p1 = i / N;
	//p2 = 1 - p1;
	//den = (sum(val(1:i)) / i - sum(val(i + 1:end)) / (N - i)) ^ 2;
	//if p1<0.5
	//	obj(i) = num / (den*p1) - p2;
	//else
	//	obj(i) = num / (den*p2) - p1;
	//end
	//	end


	//plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL); /* Create the output matrix */
	//pOut = mxGetPr(plhs[0]); /* Get the pointer to the data of B */
	//
	//Sum = mxCreateDoubleMatrix(1, N, mxREAL); 
	//pSum = mxGetPr(Sum);

	//// Calculate mean
	//for (j = 0; j != N;j++)
	//{
	//	pSum[j] = 0;
	//	for (i = 0; i != M; i++)
	//	{
	//		pSum[j]=pSum[j]+pData[i+M*j];
	//	}
	//	pSum[j]=pSum[j]/M;
	//}

	//// Subtract the mean
	//for (j = 0; j != N;j++)
	//{
	//	for (i = 0; i != M; i++)
	//	{
	//		pData[i+M*j]=pData[i+M*j]-pSum[j];
	//	}
	//}

	//// Calculate covariance
	//for (j = 0; j != N;j++)
	//{
	//	for (jj = 0; jj != N; jj++)
	//	{
	//		if (jj<j)
	//		{
	//			pOut[jj+j*N]=pOut[j+jj*N];
	//		}
	//		else
	//		{
	//			pOut[jj+j*N] = 0;
	//			for(i = 0; i != M; i++)
	//			{
	//				pOut[jj+j*N]=pOut[jj+j*N]+pData[i+M*j]*pData[i+M*jj];
	//			}
	//			pOut[jj+j*N]=pOut[jj+j*N]/M;
	//		}
	//	}
	//}

	return;
}