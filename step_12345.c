#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


int main(int argc,char *argv[])
{
	int i, j, k, nnode, *arr, flag, cnt;
	FILE *fp1, *fp2, *fp3, *fp4;
	long double **x, **sim_s_wgt, *p1_wgt, *p2_wgt, *y, **sim_s_pr, *p1_pr, *p2_pr, **sim_1, **sim_2, **w_1, **w_2, **cl_adj, sum_1_wgt, sum_1_pr, sum_1_sim, sum_2_sim, t3, t4, sum_1_cl_adj, sum_2_cl_adj, tgt_1, tgt_2, sim_wgt;
	float alpha_1, alpha_2;

	size_t len=0;
	ssize_t read;
	long int cnt_line=0;
	char delim[]=" ", *line=NULL, *ptr;



	nnode=atoi(argv[1]); //Number of nodes

	//weight matrix file name
	char *fname_1;
	fname_1=(char *)malloc(100*sizeof(char));
	const char base_1[]="./";
	const char base_2[]="/";
	sprintf(fname_1,"%s%s%s%s",base_1,argv[2],base_2,argv[3]);//Main directory name and weight matrix file name

	fp1=fopen(fname_1,"r");

	x=(long double **)malloc(nnode*sizeof(long double *));
	
	sim_s_wgt=(long double **)calloc(nnode,sizeof(long double *));

	p1_wgt=(long double *)calloc(nnode,sizeof(long double));
	p2_wgt=(long double *)calloc(nnode,sizeof(long double));

	//page rank vector file name
	char *fname_2;
	fname_2=(char *)malloc(100*sizeof(char));
	sprintf(fname_2,"%s%s%s%s",base_1,argv[2],base_2,argv[4]);//Main directory name and page rank score file name

	fp2=fopen(fname_2,"r");

	y=(long double *)malloc(nnode*sizeof(long double));

	sim_s_pr=(long double **)calloc(nnode,sizeof(long double *));

	p1_pr=(long double *)calloc(nnode,sizeof(long double));
	p2_pr=(long double *)calloc(nnode,sizeof(long double));

	sim_1=(long double **)calloc(nnode,sizeof(long double *));
	sim_2=(long double **)calloc(nnode,sizeof(long double *));


	//maximal-cliques file name
	char *fname_3;
	fname_3=(char *)malloc(100*sizeof(char));
	sprintf(fname_3,"%s%s%s%s",base_1,argv[2],base_2,argv[5]);//Main directory name and maximal-cliques file name

	fp3=fopen(fname_3,"r");


	w_1=(long double **)calloc(nnode,sizeof(long double *));
	w_2=(long double **)calloc(nnode,sizeof(long double *));


	arr=(int *)malloc(nnode*sizeof(int));


	alpha_1=atof(argv[6]);//Alpha value less than or equal to 1
	//printf("\n %f \n",alpha_1);
	alpha_2=(1-alpha_1);

	cl_adj=(long double **)calloc(nnode,sizeof(long double *));

	if(fp2==NULL)
	{
		//printf("\n Error! opening file");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(i=0;i<nnode;i++)
		{
			fscanf(fp2,"%LF\n",&y[i]);
		}
	}



	if(fp1==NULL)
	{
		//printf("\n Error! opening file");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(i=0;i<nnode;i++)
		{
			*(x+i)=(long double *)malloc(nnode*sizeof(long double));
			
			*(sim_s_wgt+i)=(long double *)calloc(nnode,sizeof(long double));

			*(sim_s_pr+i)=(long double *)calloc(nnode,sizeof(long double));


			*(sim_1+i)=(long double *)calloc(nnode,sizeof(long double));
			*(sim_2+i)=(long double *)calloc(nnode,sizeof(long double));


			*(w_1+i)=(long double *)calloc(nnode,sizeof(long double));
			*(w_2+i)=(long double *)calloc(nnode,sizeof(long double));

			*(cl_adj+i)=(long double *)calloc(nnode,sizeof(long double));

			sum_1_wgt=0;
			sum_1_pr=0;
			for(j=0;j<nnode;j++)
			{

				sim_s_wgt[i][j]=0;
				sim_s_pr[i][j]=0;

				sim_1[i][j]=0;
				sim_2[i][j]=0;

				w_1[i][j]=0;
				w_2[i][j]=0;

				if(j<=(nnode-2))
				{
					fscanf(fp1,"%LF\t",&x[i][j]);

					sum_1_wgt=sum_1_wgt+(x[i][j]*x[i][j]);

					if(x[i][j]>0)
					{
						sum_1_pr=sum_1_pr+y[j];
					}
				}
				else
				{
					fscanf(fp1,"%LF\n",&x[i][j]);

					sum_1_wgt=sum_1_wgt+(x[i][j]*x[i][j]);

					if(x[i][j]>0)
					{
						sum_1_pr=sum_1_pr+y[j];
					}
					
				}
				
			}
			p1_wgt[i]=sum_1_wgt;
			p2_wgt[i]=sum_1_wgt;

			p1_pr[i]=sum_1_pr;
			p2_pr[i]=sum_1_pr;
		}
	}


	for(i=0;i<nnode;i++)
	{
		for(j=0;j<=i;j++)
		{

			sum_1_sim=0;
			sum_2_sim=0;
	  		if(i!=j)
	  		{
				sim_wgt=sqrt(p1_wgt[i]*p2_wgt[j]);
				sim_s_wgt[i][j]=sim_wgt;
				sim_s_wgt[j][i]=sim_wgt;

				sim_s_pr[i][j]=sqrt(p1_pr[i])*sqrt(p2_pr[j]);
				sim_s_pr[j][i]=sim_s_pr[i][j];


				flag=0;
				for(k=0;k<nnode;k++)
				{
					if(x[i][k]!=0 && x[j][k]!=0)
					{
						sum_1_sim=sum_1_sim+(x[i][k]*x[j][k]);
						sum_2_sim=sum_2_sim+y[k];
						flag=1;
					}
				}
				if(flag==1)
				{
					sim_1[i][j]=x[i][j];
					sim_1[j][i]=sim_1[i][j];


					t3=sum_2_sim;
					t4=sim_s_pr[i][j];
					sim_2[i][j]=t3/t4;
					sim_2[j][i]=sim_2[i][j];
				}
				if(flag==0)
				{
					sim_1[i][j]=x[i][j];
					sim_1[j][i]=sim_1[i][j];


					t3=0;
					t4=sim_s_pr[i][j];
					sim_2[i][j]=t3/t4;
					sim_2[j][i]=sim_2[i][j];
				}
				
	  		}
		}
	}

	
	if(fp3==NULL)
	{
		//printf("\n Error! opening file");
		exit(EXIT_FAILURE);
	}
	else
	{
		while((read=getline(&line,&len,fp3))!=-1)
		{
			cnt_line++;
			if(cnt_line>2)
			{
				ptr=strtok(line,delim);
				cnt=0;
				while(ptr!=NULL)
				{
					//printf("%d\n",atoi(ptr));
					arr[cnt]=atoi(ptr);
					cnt++;
					ptr=strtok(NULL,delim);
				}

				sum_1_cl_adj=0;
				sum_2_cl_adj=0;

        			for(j=0;j<cnt;j++)
        			{
            				//sum_3=0;
            				for(k=0;k<cnt;k++)
            				{
                				if(j!=k)
                				{
                    					sum_1_cl_adj=sum_1_cl_adj+sim_1[arr[j]][arr[k]];
                    					sum_2_cl_adj=sum_2_cl_adj+sim_2[arr[j]][arr[k]];
                				}
           		 		}
        			}
				tgt_1=(sum_1_cl_adj/(cnt*(cnt-1)));
				tgt_2=(sum_2_cl_adj/(cnt*(cnt-1)));

				for(j=0;j<cnt;j++)
				{
					for(k=0;k<cnt;k++)
					{
						if(j!=k)
						{
							w_1[arr[j]][arr[k]]=w_1[arr[j]][arr[k]]+tgt_1;
							w_2[arr[j]][arr[k]]=w_2[arr[j]][arr[k]]+tgt_2;
						}
					}
				}

			}
		}
	}


	char *fname_4;
	fname_4=(char *)malloc(100*sizeof(char));
	sprintf(fname_4,"%s%s%s%s",base_1,argv[2],base_2,argv[7]);//Main directory name and clique-adjacency matrix file name


	fp4=fopen(fname_4,"w");

	if(fp4==NULL)
	{
		//printf("\n Error! opening file");
		exit(EXIT_FAILURE);
	}
	else
	{
		for(i=0;i<nnode;i++)
		{
			for(j=0;j<nnode;j++)
			{
				if(j<=(nnode-2))
				{
					cl_adj[i][j]=(alpha_1*w_1[i][j])+(alpha_2*w_2[i][j]);
					fprintf(fp4,"%LF\t",cl_adj[i][j]);
				}
				else
				{
					cl_adj[i][j]=(alpha_1*w_1[i][j])+(alpha_2*w_2[i][j]);
					fprintf(fp4,"%LF\n",cl_adj[i][j]);
				}
			}
		}
	}


	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	free(x);
	free(sim_s_wgt);
	free(p1_wgt);
	free(p2_wgt);
	free(y);
	free(sim_s_pr);
	free(p1_pr);
	free(p2_pr);
	free(sim_1);
	free(sim_2);
	free(w_1);
	free(w_2);
	free(cl_adj);


	free(fname_1);
	free(fname_2);
	free(fname_3);
	free(fname_4);
		

	return 0;
}
				
