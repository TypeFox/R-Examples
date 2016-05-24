#
#  Copyright (C) 2007 Sebastian Kaiser
#  Bicluster Algorithm for Questionairs based on Murali, T. & Kasif, S. Extracting conserved gene expression motifs from gene expression data Proc. Pacific Symp. Biocomputing, sullivan.bu.edu, 2003
#

## algorithm to find the biggest bicluster in questionaire (modified xmotif algorithm for questionairs)

bigquestordmotif<-function(mat,d,ns,nd,sd,alpha)
{
    if(length(d)==1) d[2] <- d
size<-4
nr<-nrow(mat)
ques<- rep(FALSE,ncol(mat))
pers<- rep(FALSE,nrow(mat))
for(i in 1:ns)
{
ri<-sample(1:nr,1)
logr<-rep(TRUE,nrow(mat))
logr[ri]<-FALSE

for(j in 1:nd)
{
D<-sample(1:nr,sd,prob=logr)

gri<-mat[ri,]
griD<-c(D,ri)
cS<-rowSums((t(mat[griD,])>=(gri-d[1])) & (t(mat[griD,])<=(gri+d[2])))
gij<-cS==length(griD)

if(sum(gij)>= max(sum(ques),2))
{
rri<-mat[ri,gij]
rS<-colSums((t(mat[,gij])>=(rri-d[1])) & (t(mat[,gij])<=(rri+d[2])))
rij<-rS==sum(gij)

if((sum(rij)>=(alpha*nr)) & ((sum(gij)*sum(rij))>size) )
{
ques<-gij
pers<-rij
size <- sum(ques)*sum(pers)
}
}
}

}
erg<-list(pers,ques)
erg
}




## algorithm to find number biggest bicluster (Stops if all persons are in one bicluster or if no bicluster is found)

questordmotif<-function(mat,d,ns,nd,sd,alpha,number)
{
MYCALL <- match.call()
x<-matrix(FALSE,nrow=nrow(mat),ncol=number)
y<-matrix(FALSE,nrow=number,ncol=ncol(mat))
matstore<-mat
Stop <- FALSE
logr<-rep(TRUE,nrow(mat))
for(i in 1:number)
{
erg<-bigquestordmotif(mat,d,ns,nd,sd,alpha)
if(sum(erg[[1]])==0)
{
    Stop <- TRUE
    break
}
else{
x[logr,i]<-erg[[1]]
y[i,]<-erg[[2]]
logr[logr][erg[[1]]]<-FALSE
mat<-matstore[logr,]
if(sum(logr)<(sd+1))
{
    Stop <- TRUE
    break}
}
}
if(Stop)
{return(BiclustResult(as.list(MYCALL),as.matrix(x[,1:(i-1)]),as.matrix(y[1:(i-1),]),(i-1),list(0)))
}
else{
return(BiclustResult(as.list(MYCALL),as.matrix(x),as.matrix(y),i,list(0)))
}
}
