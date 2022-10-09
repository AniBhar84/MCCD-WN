#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>


void normalize(long double **cl_data,long int nnode,long int n_cm,long double **norm_cl_data,long double *rnorm);
void dist_func(long double **norm_cl_data, long int nnode,long int n_cm,long double **man_dist,long double **eucl_dist);
void center_select(long int nnode,long int n_cm,long double **eucl_dist,long int *ind,long int sid);
void write(long int *clustering,long int nnode, char *argv[]);
long double A_kval_comp(long int **x,long double *wgt, long int nnode,long double **A,long double *kval, long int count);
long double mod_comp(long double **A,long double *kval,long int n_cm,long int nnode,long double wgt_sum,long int *clustering);
void ind_mat_gen(long double **norm_cl_data,long int n_cm,long int *ind,long double **ind_mat);
void dist_mat_func(long double **norm_cl_data,long int nnode, long int n_cm,long double **ind_mat,long double **eu_dist_mat,long double **man_dist_mat);
void kmed(long double **norm_cl_data,long double **eu_dist_mat,long double **man_dist_mat,long int nnode,long int n_cm,long int *clustering,long double *bwsc,long double **ind_mat);
void wz_algo(long double **norm_cl_data,long int nnode,long int n_cm,long double **eu_dist_mat,long double **man_dist_mat,long int *clustering,long double **ind_mat);
long double l2norm(long double *Ru,long int n_cm);
long double dense_score(long double **A,long double *kval,long int n_cm,long int nnode,long int *clustering);
int main(int argc,char *argv[])
{
	long int i, j, nnode, sid, *ind, K, n_cm, min_ind, *clustering, l, m, **x, Qmax_ind, K1, K2;
	long int count, i1, j1;
	long double **lap_mat, **bcup_mat, **cl_data, *bwsc, min, **man_dist, **eucl_dist, **norm_cl_data, *rnorm, sum, *wgt, **A, *kval, wgt_sum, Q, Qmax, **ind_mat, **eu_dist_mat, **man_dist_mat;
	FILE *fp1, *fp2, *fp3;
	char c;

	nnode=atol(argv[1]);
	K=atol(argv[2]);

	lap_mat=(long double **)malloc(nnode*sizeof(long double *));
	bcup_mat=(long double **)calloc(nnode,sizeof(long double *));
	bwsc=(long double *)calloc(nnode,sizeof(long double));
	kval=(long double *)calloc(nnode,sizeof(long double));

	
	

	eucl_dist=(long double **)calloc(nnode,sizeof(long double *));
	man_dist=(long double **)calloc(nnode,sizeof(long double *));
	A=(long double **)calloc(nnode,sizeof(long double *));

	char *fname_1;
	fname_1=(char *)malloc(100*sizeof(char));
	const char base_1[]="./";
	
	const char base_2[]="/";
	sprintf(fname_1,"%s%s%s%s",base_1,argv[5],base_2,argv[6]);//Main directory name and laplacian matrix file name

	char *fname_2;
	fname_2=(char *)malloc(100*sizeof(char));
	sprintf(fname_2,"%s%s%s%s",base_1,argv[5],base_2,argv[7]);//Main directory name and betweenness score file name

	
	fp1=fopen(fname_1,"r");
	fp2=fopen(fname_2,"r");
	
	if(fp1==NULL || fp2==NULL)
	{
		printf("\n File opening error 1!!");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(i=0;i<nnode;i++)
		{
			*(lap_mat+i)=(long double *)malloc(K*sizeof(long double));
			*(bcup_mat+i)=(long double *)calloc(K,sizeof(long double));
			

			*(eucl_dist+i)=(long double *)malloc(nnode*sizeof(long double));
			*(man_dist+i)=(long double *)calloc(nnode,sizeof(long double));
			*(A+i)=(long double *)calloc(nnode,sizeof(long double));
			fscanf(fp2,"%LF\n",&bwsc[i]);
			if(i==0)
			{
				min=bwsc[i];
				min_ind=i;
			}
			else
			{
				if(bwsc[i]<min)
				{
					min=bwsc[i];
					min_ind=i;
				}
			}
			for(j=0;j<K;j++)
			{
				if(j<=(K-2))
				{
					fscanf(fp1,"%LF\t",&lap_mat[i][j]);
					/*if(i==0)
						printf("%.15LF\t",lap_mat[i][j]);*/
					bcup_mat[i][j]=lap_mat[i][j];
				}
				else
				{
					fscanf(fp1,"%LF\n",&lap_mat[i][j]);
					/*if(i==0)
						printf("%.15LF\t",lap_mat[i][j]);*/
					bcup_mat[i][j]=lap_mat[i][j];	
				}
			}
		}
	}
	sid=min_ind;
	
	

	
	char *fname_3;
	fname_3=(char *)malloc(100*sizeof(char));
	sprintf(fname_3,"%s%s",base_1,argv[8]);//Main directory name and Weighted Network Edge Table file name

	fp3=fopen(fname_3,"r");

	count=0;
	if(fp3==NULL)
	{
		printf("\n File opening error 2!!");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(c=getc(fp3);c!=EOF;c=getc(fp3))
		{
			if(c=='\n')
			{
            			count=count+1;
			}
		}

	}
	
	fclose(fp3);

	
	char *fname_4;
	fname_4=(char *)malloc(100*sizeof(char));
	
	sprintf(fname_4,"%s%s",base_1,argv[8]);//Main directory name and Weighted Network Edge Table file name

	
	fp3=fopen(fname_4,"r");

	x=(long int **)malloc(count*sizeof(long int *));
	wgt=(long double *)malloc(count*sizeof(long double));

	if(fp3==NULL)
	{
		printf("\n File opening error 3!!");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(i1=0;i1<count;i1++)
		{
			*(x+i1)=(long int *)malloc(2*sizeof(long int));
			for(j1=0;j1<3;j1++)
			{
				if(j1<=(3-2))
				{
					fscanf(fp3,"%ld\t",&x[i1][j1]);
				}
				else
				{
					fscanf(fp3,"%LF\n",&wgt[i1]);	
				}
			}
		}	

	}
	fclose(fp3);
	printf("\n Reading Done!!!");
	
	wgt_sum=A_kval_comp(x,wgt,nnode,A,kval,count);

	rnorm=(long double *)calloc(nnode,sizeof(long double));
	clustering=(long int *)calloc(nnode,sizeof(long int));
	
	K1=atoi(argv[3]);
	K2=atoi(argv[4]);
	for(n_cm=K1;n_cm<=K2;n_cm++)
	{
		cl_data=(long double **)malloc(nnode*sizeof(long double *));
		norm_cl_data=(long double **)malloc(nnode*sizeof(long double *));

		
		if(n_cm==2)
		{
			
			ind=(long int *)calloc(n_cm,sizeof(long int));
			for(i=0;i<nnode;i++)
			{
				*(cl_data+i)=(long double *)malloc((n_cm-1)*sizeof(long double));
				*(norm_cl_data+i)=(long double *)malloc((n_cm-1)*sizeof(long double));
				for(j=K-n_cm,l=0;j>=K-n_cm;j--,l++)
				{
					cl_data[i][l]=lap_mat[i][j];
					
					if(cl_data[i][l]>0)
					{
						clustering[i]=1;
					}
					else
					{
						clustering[i]=2;
					}
				}
			}
			Q=mod_comp(A,kval,n_cm,nnode,wgt_sum,clustering);
			
			Qmax=Q;
			Qmax_ind=n_cm;
			
			write(clustering,nnode,argv);
			free(cl_data);
			free(norm_cl_data);
		}
		else
		{
			
			ind=(long int *)malloc(n_cm*sizeof(long int));
			ind_mat=(long double **)malloc(n_cm*sizeof(long double *));
			eu_dist_mat=(long double **)calloc(nnode,sizeof(long double *));
			man_dist_mat=(long double **)calloc(nnode,sizeof(long double *));
			for(i=0;i<nnode;i++)
			{
				*(cl_data+i)=(long double *)malloc(n_cm*sizeof(long double));
				*(norm_cl_data+i)=(long double *)malloc(n_cm*sizeof(long double));
				*(eu_dist_mat+i)=(long double *)malloc(n_cm*sizeof(long double));
				*(man_dist_mat+i)=(long double *)malloc(n_cm*sizeof(long double));
				sum=0;
				for(j=K-1,l=0;j>=K-n_cm;j--,l++)
				{
					*(ind_mat+l)=(long double *)malloc(n_cm*sizeof(long double));
					cl_data[i][l]=lap_mat[i][j];
					sum=sum+(cl_data[i][l]*cl_data[i][l]);
				}
				rnorm[i]=sqrt(sum);
				//printf("%LF \t",rnorm[i]);
			}
			
			
			normalize(cl_data,nnode,n_cm,norm_cl_data,rnorm);
			dist_func(norm_cl_data,nnode,n_cm,man_dist,eucl_dist);
			center_select(nnode,n_cm,eucl_dist,ind,sid);
			ind_mat_gen(norm_cl_data,n_cm,ind,ind_mat);
			
			
			dist_mat_func(norm_cl_data,nnode,n_cm,ind_mat,eu_dist_mat,man_dist_mat);
			kmed(norm_cl_data,eu_dist_mat,man_dist_mat,nnode,n_cm,clustering,bwsc,ind_mat);
			Q=mod_comp(A,kval,n_cm,nnode,wgt_sum,clustering);
			
			if(Q>Qmax)
			{
				Qmax=Q;
				Qmax_ind=n_cm;
			}
			write(clustering,nnode,argv);
			free(cl_data);
			free(norm_cl_data);
			free(eu_dist_mat);
			free(man_dist_mat);
			free(ind_mat);
		}
		
	}
	if(K1!=K2)
	{
		printf("\n The no. of communities is %ld and max. modularity is %LF \n",Qmax_ind,Qmax);
	}
	fclose(fp1);
	fclose(fp2);
	
	free(ind);
	free(clustering);
	free(lap_mat);
	free(bcup_mat);
	free(bwsc);
	free(man_dist);
	free(eucl_dist);
	free(rnorm);
	free(x);
	free(A);
	free(wgt);
	free(kval);

	free(fname_1);
	free(fname_2);
	free(fname_3);
	free(fname_4);
	
	return(0);
}
		
void ind_mat_gen(long double **norm_cl_data,long int n_cm,long int *ind,long double **ind_mat)
{
	long int i, j;
	for(i=0;i<n_cm;i++)
	{
		for(j=0;j<n_cm;j++)
		{
			ind_mat[i][j]=norm_cl_data[ind[i]][j];
		}
	}


}


void dist_mat_func(long double **norm_cl_data,long int nnode, long int n_cm,long double **ind_mat,long double **eu_dist_mat,long double **man_dist_mat)
{
	long int i, j, k;
	long double ds_eu, ds_man, d;
	for(i=0;i<nnode;i++)
	{
		for(j=0;j<n_cm;j++)
		{
			ds_eu=0;
			ds_man=0;
			for(k=0;k<n_cm;k++)
			{
				d=(norm_cl_data[i][k]-ind_mat[j][k]);
				if(d<0)
				{
					d=(-1.0)*d;
					ds_eu=ds_eu+(d*d);
					ds_man=ds_man+d;
				}
				else
				{
					d=(1.0)*d;
					ds_eu=ds_eu+(d*d);
					ds_man=ds_man+d;
				}
			}
			eu_dist_mat[i][j]=sqrt(ds_eu);
			man_dist_mat[i][j]=ds_man;
		}
	}	
}

void normalize(long double **cl_data,long int nnode,long int n_cm,long double **norm_cl_data,long double *rnorm)
{
	long int i, j;
	for(i=0;i<nnode;i++)
	{	
		for(j=0;j<n_cm;j++)
		{
			norm_cl_data[i][j]=cl_data[i][j]/rnorm[i];

		} 
		

	}

}

void dist_func(long double **norm_cl_data, long int nnode,long int n_cm,long double **man_dist,long double **eucl_dist)
{
	long int i, j, k;
	long double s0, s1, s2;
	for(i=0;i<nnode;i++)
	{
		for(j=0;j<=i;j++)
		{
			s1=0;
			s2=0;
			if(i==j)
			{
				man_dist[i][j]=0;
				eucl_dist[i][j]=0;
				man_dist[j][i]=0;
				eucl_dist[j][i]=0;
			}
			else
			{
				for(k=0;k<n_cm;k++)
				{
					s0=(norm_cl_data[i][k]-norm_cl_data[j][k]);
					if(s0<0)
					{
						s1=s1+((-1.0)*s0);
					}
					if(s0>=0)
					{
						s1=s1+((1.0)*s0);
					}
					s2=s2+((norm_cl_data[i][k]-norm_cl_data[j][k])*(norm_cl_data[i][k]-norm_cl_data[j][k]));
				}
				man_dist[i][j]=s1;
				eucl_dist[i][j]=sqrt(s2);
				man_dist[j][i]=s1;
				eucl_dist[j][i]=sqrt(s2);
			}
		}
	}

}


void center_select(long int nnode,long int n_cm,long double **eucl_dist,long int *ind,long int sid)
{
	long int t=1, i, j, flag, min_ind, cnt1, max_ind;
	long double min, max;
	ind[0]=sid;
	while(t<n_cm)
	{
		cnt1=0;
		for(i=0;i<nnode;i++)
		{
			flag=1;
			for(j=0;j<t;j++)
			{
				if(ind[j]==i)
				{
					flag=0;
					break;
					
				}
			}
			if(flag==0)
			{
				continue;
			}
			if(flag==1)
			{
				for(j=0;j<t;j++)
				{
					if(j==0)
					{
						min=eucl_dist[i][ind[j]];
						min_ind=ind[j];
					}
					else
					{
						if(eucl_dist[i][ind[j]]<min)
						{
							min=eucl_dist[i][ind[j]];
							min_ind=ind[j];
						}
					}
				}
				if(cnt1==0)
				{
					max=(min*min);
					max_ind=i;
					cnt1++;
				}
				if(cnt1>0)
				{
					if((min*min)>max)
					{
						max=(min*min);
						max_ind=i;
						cnt1++;
					}
				}
				
			}
		}
		ind[t]=max_ind;
		t++;
		
	}

}


void kmed(long double **norm_cl_data,long double **eu_dist_mat,long double **man_dist_mat,long int nnode,long int n_cm,long int *clustering,long double *bwsc,long double **ind_mat)
{
	long int i, j, k, min_pos,iter=0, flag=1, fl, *cl_bcup;
	long double min_ds;
	cl_bcup=(long int *)calloc(nnode,sizeof(long int));
	while(flag)
	{
		fl=0;
		iter++;
		for(i=0;i<nnode;i++)
		{
			for(k=0;k<n_cm;k++)
			{
				if(k==0)
				{
					min_ds=eu_dist_mat[i][k];
					min_pos=k;
				}
				else
				{
					if(eu_dist_mat[i][k]<min_ds)
					{
						min_ds=eu_dist_mat[i][k];
						min_pos=k;
					}
				}
			}
			clustering[i]=(min_pos+1);
			if(iter==1)
			{
				cl_bcup[i]=clustering[i];
				fl=1;
			}
			else
			{
				if(cl_bcup[i]!=clustering[i])
				{
					fl=1;
				}
				cl_bcup[i]=clustering[i];
			}
		}
		if(fl==1)
		{
			wz_algo(norm_cl_data,nnode,n_cm,eu_dist_mat,man_dist_mat,clustering,ind_mat);
			dist_mat_func(norm_cl_data,nnode,n_cm,ind_mat,eu_dist_mat,man_dist_mat);
			flag=1;
		}
		if(fl==0)
		{
			flag=0;
		}
		
	}
	free(cl_bcup);
	
}


void wz_algo(long double **norm_cl_data,long int nnode,long int n_cm,long double **eu_dist_mat,long double **man_dist_mat,long int *clustering,long double **ind_mat)
{

	long int i, j, k;
	long double *num, wgt, wgt_sum, etau, *Ru, runorm, gamau;
	
	for(k=0;k<n_cm;k++)
	{
		
		wgt_sum=0;
		num=(long double *)calloc(n_cm,sizeof(long double));
		etau=0;
		Ru=(long double *)calloc(n_cm,sizeof(long double));
		for(i=0;i<nnode;i++)
		{
			if(clustering[i]==(k+1))
			{
				if(eu_dist_mat[i][k]!=0)
				{
					//wgt=1;
					//wgt_sum=wgt_sum+wgt;
					wgt=1/eu_dist_mat[i][k];
					wgt_sum=wgt_sum+wgt;
					for(j=0;j<n_cm;j++)
					{
						num[j]=num[j]+(norm_cl_data[i][j]*wgt);
						Ru[j]=Ru[j]+((norm_cl_data[i][j]-ind_mat[k][j])*wgt);
					}
				}
				if(eu_dist_mat[i][k]==0)
				{
					etau=1;
				}
				
			
			}
		}
		runorm=l2norm(Ru,n_cm);
		if(1<(etau/runorm))
		{
			gamau=1;
		}
		else
		{
			gamau=(etau/runorm);
		}
		
		for(j=0;j<n_cm;j++)
		{
			ind_mat[k][j]=((1-gamau)*(num[j]/wgt_sum))+((gamau)*ind_mat[k][j]);
		}
		free(num);
	}



}

long double l2norm(long double *Ru,long int n_cm)
{
	long int i;
	long double norm=0;
	for(i=0;i<n_cm;i++)
	{
		norm=norm+(Ru[i]*Ru[i]);
	}

	return(sqrt(norm));

}

void write(long int *clustering,long int nnode, char *argv[])
{
	FILE *fp4;
	long int i;

	char *fname_5;
	fname_5=(char *)malloc(100*sizeof(char));
	const char base_3[]="./";
	const char base_4[]="/comm.txt";
	sprintf(fname_5,"%s%s%s",base_3,argv[5],base_4);

	fp4=fopen(fname_5,"w");

	for(i=0;i<nnode;i++)
	{
		fprintf(fp4,"%ld\n",clustering[i]);
	
	}
	fclose(fp4);
	free(fname_5);

}

long double A_kval_comp(long int **x,long double *wgt, long int nnode,long double **A,long double *kval, long int count)
{
	long int i, tar;
	long int j;
	long double s1, s2=0;
	for(i=0;i<nnode;i++)
	{
		s1=0;
		for(j=0;j<count;j++)
		{
			if(x[j][0]==(i+1))
			{
				tar=x[j][1];
				A[i][tar-1]=wgt[j];
				s1=s1+wgt[j];
			}
		}
		kval[i]=s1;
		s2=s2+s1;
	}
	return(s2);
}


long double mod_comp(long double **A,long double *kval,long int n_cm,long int nnode,long double wgt_sum,long int *clustering)
{
	long int i, j;
	long double mvv, mval=0, Q;
	for(i=0;i<nnode;i++)
	{
		for(j=0;j<=i;j++)
		{
			if(i!=j)
			{
				if(clustering[i]==clustering[j])
				{
					mvv=(A[i][j]-((kval[i]*kval[j])/(wgt_sum)));
					mval=mval+(mvv+mvv);
				}
			}
		}
	}
	Q=(1/(wgt_sum))*mval;
	return(Q);
}


long int pam_swap_new(long double **man_dist,long int n_cm,long int nnode,long int *ind,long int *clustering,long double *bwsc)
{
	long int i, k, f, m_b_pos, flag=0;
	long double m_b;
	for(k=0;k<n_cm;k++)
	{
		f=0;
		for(i=0;i<nnode;i++)
		{
			if(clustering[i]==(k+1))
			{
				if(f==0)
				{
					m_b=bwsc[i];
					m_b_pos=i;
					f++;
				}
				else
				{
					if(bwsc[i]<m_b)
					{
						m_b=bwsc[i];
						m_b_pos=i;
						f++;

					}
				}
			}
		}
		if(ind[k]!=m_b_pos)
		{
			ind[k]=m_b_pos;
			flag=1;
		}

	}
	if(flag==1)
	{
		return(flag);
	}
	if(flag==0)
	{
		return(flag);
	}

}


long double dense_score(long double **A,long double *kval,long int n_cm,long int nnode,long int *clustering)
{
	long int i, j, k;
	long double deno, num, ds=0;
	for(k=1;k<=n_cm;k++)
	{
		
		for(i=0;i<nnode;i++)
		{
			if(clustering[i]==k)
			{
				
				deno=kval[i];
				num=0;
				for(j=0;j<nnode;j++)
				{
					if(clustering[j]==k)
					{
						num=num+A[i][j];
					}
				}
				ds=ds+(num/deno);
			}
		}
	}
	return(ds);
}
					
