#
#  Copyright (C) 2007 Sebastian Kaiser
#  Murali, T. & Kasif, S. Extracting conserved gene expression motifs from gene expression data Proc. Pacific Symp. Biocomputing, sullivan.bu.edu, 2003
#

## algorithm to find the biggest bicluster (directly programmed FindMotif algorithm of the paper)

bigxmotif<-function(mat,ns,nd,sd,alpha)
{
#Preprocess Data step not yet implemented
#statemat<-makestatemat(mat)
size<-4
nc<-ncol(mat)
gen<- rep(FALSE,nrow(mat))
co<- rep(FALSE,ncol(mat))
for(i in 1:ns)
{
ci<-sample(1:nc,1)
logc<-rep(TRUE,ncol(mat))
logc[ci]<-FALSE

for(j in 1:nd)
{
D<-sample(1:nc,sd,prob=logc)

gci<-mat[,ci]
gciD<-c(D,ci)
rS<-rowSums(mat[,gciD]==gci)
gij<-rS==length(gciD)

if(sum(gij) >= max(sum(gen),2))
{
cci<-mat[gij,ci]
cS<-colSums(mat[gij,]==cci)
cij<-cS==sum(gij)

if(sum(cij)>=(alpha*nc) & ((sum(gij)*sum(cij))>size) )
{
gen<-gij
co<-cij
size <- sum(gij) * sum(cij)
}

}
}

}
erg<-list(gen,co)
erg
}


## algorithm which stores all find motifs not only the biggest

storexmotif<-function(mat,ns,nd,sd,alpha)
{
#Preprocess Data step not yet implemented
#statemat<-makestatemat(mat)
nc<-ncol(mat)
xstore<-matrix(0,nrow=nrow(mat),ncol=(ns*nd))
ystore<-matrix(0,nrow=(ns*nd),ncol=ncol(mat))
for(i in 1:ns)
{
ci<-sample(1:nc,1)
logc<-rep(TRUE,ncol(mat))
logc[ci]<-FALSE

for(j in 1:nd)
{
D<-sample(1:nc,sd,prob=logc)

gci<-mat[,ci]
gciD<-c(D,ci)
rS<-rowSums(mat[,gciD]==gci)
gij<-rS==length(gciD)
if(sum(gij)>2)
{

cci<-mat[gij,ci]
cS<-colSums(mat[gij,]==cci)
cij<-cS==sum(gij)

xstore[,j+(i-1)*nd]<-gij*1
ystore[j+(i-1)*nd,]<-cij*1
}
}

}
erg<-list(xstore,ystore)
erg


}


## algorithm to find number biggest bicluster (Stops if all genes are in one bicluster or if no bicluster is found)

xmotifbiclust<-function(mat,ns,nd,sd,alpha,number)
{
MYCALL <- match.call()
x<-matrix(FALSE,nrow=nrow(mat),ncol=number)
y<-matrix(FALSE,nrow=number,ncol=ncol(mat))
matstore<-mat
Stop <- FALSE
logr<-rep(TRUE,nrow(mat))
for(i in 1:number)
{
    if(sum(logr) < 2)
        {
            Stop <- TRUE
            break
        }
    erg<-bigxmotif(mat,ns,nd,sd,alpha)
    if(sum(erg[[1]])==0)
    {
        Stop <- TRUE
        break
    }
    else
    {
        x[logr,i]<-erg[[1]]
        y[i,]<-erg[[2]]
        logr[logr][erg[[1]]]<-FALSE
        mat<-as.matrix(matstore[logr,])
    }
}
if(Stop)
{
    return(BiclustResult(as.list(MYCALL),as.matrix(x[,1:(i-1)]),as.matrix(y[1:(i-1),]),(i-1),list(0)))
}
else
{
    return(BiclustResult(as.list(MYCALL),as.matrix(x),as.matrix(y),i,list(0)))
}
}
