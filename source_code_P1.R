#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
tp<-as.character(args[1])#Main directory for output files

library(data.table)
library(igraph)

#fname_1<-"./novl_0.nsa"
fname_1<-as.character(args[2])#Weighted Network Edge Table file name (tab delimited)
x_wgt<-fread(fname_1,sep="\t",header=FALSE)
x_wgt<-as.matrix(x_wgt)

node_name<-c(1:length(unique(as.numeric(x_wgt[,1]))))

#fname_2<-paste("./",tp,"/node_name.txt",sep="")
#write.table(node_name,fname_2,append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)

x_wgt_data<-matrix(0,length(node_name),length(node_name))
for(i in 1:length(node_name))
{
	pos_i<-which(as.numeric(x_wgt[,1])==i)
	tar_i<-as.numeric(x_wgt[pos_i,2])
	wgt_i<-as.numeric(x_wgt[pos_i,3])
	x_wgt_data[i,tar_i]<-wgt_i
}

st_1<-as.character(args[3])#Output file name for weighted adjacency matrix
#fname_3<-paste("./",tp,"/ori_weight.txt",sep="")
fname_3<-paste("./",tp,"/",st_1,sep="")
write.table(x_wgt_data,fname_3,sep="\t",append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)

if(FALSE)#for simulated network
{
	ind_nw_final<-NULL
	for(i in 1:nrow(x_wgt_data))
	{
		tar_ind<-which(as.numeric(x_wgt_data[i,])>0)
		c1<-matrix(rep(i,times=length(tar_ind)),length(tar_ind),1)	
		c2<-matrix(tar_ind,length(tar_ind),1)
		c3<-matrix(as.numeric(x_wgt_data[i,tar_ind]),length(tar_ind),1)
		ind_nw_final<-rbind(ind_nw_final,cbind(cbind(c1,c2),c3))

	}
	fname_4<-paste("./",tp,"/ind_nw_final.txt",sep="")
	write.table(ind_nw_final,fname_4,sep="\t",append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)
}
#*************Part 2***********
n_edges<-0
nw_final<-NULL
for(i in 1:nrow(x_wgt_data))
{
	pos_nzero<-which(as.numeric(x_wgt_data[i,])>0)
	c1<-rep(i-1,times=length(pos_nzero))
	c1<-matrix(c1,length(c1),1)
	c2<-(pos_nzero-1)
	c2<-matrix(c2,length(c2),1)
	nw_final<-rbind(nw_final,cbind(c1,c2))
	n_edges<-n_edges+length(pos_nzero)
}
n_vertex<-nrow(x_wgt_data)
s1<-append(n_vertex,n_edges)

st_2<-as.character(args[4])#Output file name for edge list (maximal-clique algorithm)
#fname_6<-paste("./",tp,"/edge_list.txt",sep="")
fname_6<-paste("./",tp,"/",st_2,sep="")
write.table(s1,fname_6,append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(nw_final,fname_6,sep=",",append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)




x_wgt<-as.data.frame(x_wgt)
g_1<-graph_from_edgelist(as.matrix(x_wgt[,c(1,2)]),directed=TRUE)
E(g_1)$weight<-as.numeric(x_wgt[,3])

ks<-page_rank(g_1,algo="prpack",vids=V(g_1),directed=F,damping=0.85,weights=1/E(g_1)$weight,options=NULL)

st_3<-as.character(args[5])#Output file name for page-rank score of vertices
#fname_7<-paste("./",tp,"/pr_score.txt",sep="")
fname_7<-paste("./",tp,"/",st_3,sep="")
write.table(as.numeric(ks$vector),fname_7,append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)


