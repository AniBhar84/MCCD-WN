#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
tp<-as.character(args[1])#Main directory for output files

MAX_K<-as.numeric(args[2])#Maximum number of eigen-vectors to be considered

suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

suppressPackageStartupMessages(library(rARPACK))


st_1<-as.character(args[3])#Clique-adjacency matrix file name
fname_2<-paste("./",tp,"/",st_1,sep="")
w_1<-read.delim(fname_2,sep="\t",header=FALSE)
w_1<-as.data.frame(w_1)
w_1<-apply(w_1,2,as.numeric)
w_1<-as.matrix(w_1)


#*****weights of nodes*****
suppressPackageStartupMessages(library(igraph))
st_2<-as.character(args[4])#Weighted Network Edge Table file name (tab delimited)
fname_3<-st_2
x<-read.delim(fname_3,sep="\t",header=FALSE)
#x<-as.data.frame(x)
#x<-apply(x,2,as.numeric)
#x<-as.matrix(x)

g_1<-graph_from_edgelist(as.matrix(x[,c(1,2)]),directed=TRUE)
E(g_1)$weight<-as.numeric(x[,3])
bwsc<-betweenness(g_1,v=V(g_1),directed=TRUE,weights=1/E(g_1)$weight,normalized=T)

st_3<-as.character(args[5])#Betweenness centrality score of vertices
fname_4<-paste("./",tp,"/",st_3,sep="")
write.table(as.data.frame(bwsc),fname_4,append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)



#************************Part 3************************
deg_1<-as.numeric(rowSums(w_1))
half_d<-diag(deg_1^(-1/2))
deg_1<-diag(deg_1)


l_1<-(half_d%*%(deg_1-w_1)%*%half_d)



max_col<-MAX_K
eig_lap<-eigs(l_1,max_col,which="SM",sigma=NULL)
data_cl<-eig_lap$vectors
bup_cl<-data_cl
data_cl_val<-eig_lap$values

st_4<-as.character(args[6])#output laplacian matrix file name

fname_5<-paste("./",tp,"/",st_4,sep="")
write.table(data_cl,fname_5,sep="\t",append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)



