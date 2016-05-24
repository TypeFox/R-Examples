# createAdj.R
# create an adjancency matrix by first creating a matrix

createAdj<-function(matrix,filename='Adj.txt',suffix=NULL){
# checks
if(is.matrix(matrix)!=TRUE){stop('Input must be a matrix')}
if(dim(matrix)[1]!=dim(matrix)[2]){stop('Matrix must be square')}
if(isSymmetric(matrix)==FALSE){stop('Matrix must be symmetric')}
# Vectors from matrix
n<-nrow(matrix)
num<-vector(length=n,mode='numeric')
for (i in 1:n){
   num[i]<-sum(matrix[i,],na.rm=TRUE)
}
adj<-vector(length=sum(num),mode='numeric')
weight<-vector(length=sum(num),mode='numeric')
index<-1
for (i in 1:n){
   ind<- is.na(matrix[i,])==FALSE
   xxx<-as.numeric(ind)*(1:n)
   if (sum(xxx)>0){
      aaa<-xxx[ind]
      www<-matrix[i,ind]
      adj[index:(index+length(aaa)-1)]<-aaa
      weight[index:(index+length(aaa)-1)]<-www
      index<-index+length(aaa)
   }
}
sumNumNeigh<-sum(weight)
## Note for consideration next version
## consider using something like
## dump(ls(pattern = paste("num",suffix,sep="")),file=filename)
## although no need to change if this is working OK     PB 4/12/2009
## Create adjacency matrix in CAR format ##
zz <- file(filename, "w")  # open an output file connection
cat("list(num",suffix,"=c(\n",sep='', file = zz)
nums<-paste(num,collapse=',')
cat(nums,sep='',file = zz)
cat("),\nadj",suffix,"=c(\n",sep='', file = zz)
# Output adjacency numbers as one row per region
index<-0
for (i in 1:length(num)){
   if (num[i]>0){
     for (j in 1:num[i]){
        index=index+1;
        if(index<sumNumNeigh){cat(adj[index],",",sep='',file = zz)}
        if(index==sumNumNeigh){cat(adj[index],sep='',file = zz)}
     }
   }
   cat("\n",sep='',file = zz)
}
cat("),\nweights",suffix,"=c(\n",sep='', file = zz)
weights<-paste(weight,collapse=',')
cat(weights,sep='',file = zz)
cat("))\n",sep='', file = zz)
close(zz)
toret<-list()
toret$num<-num
toret$adj<-adj
toret$weight<-weight
return(toret)
}
# Example (nearest neighbour)
#x<-c(NA,1,NA,NA,NA)
#V= toeplitz (x)
#createAdj(V)
