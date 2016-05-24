#
#  Copyright (C) 2007 Sebastian Kaiser
#  Cheng, Y. & Church, G. Biclustering of expression data. Proc Int Conf Intell Syst Mol Biol, ncbi.nlm.nih.gov, 2000
#


##Some helper funktions to calculate the CC score for node deletion

ccscore<-function(mat)
{
score<-sum((mat-rowMeans(mat)-matrix(colMeans(mat),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat))^2)/(nrow(mat)*ncol(mat))
score
}

rowscore<-function(mat)
{
score<-rowSums((mat-rowMeans(mat)-matrix(colMeans(mat),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat))^2)/ncol(mat)
score
}

colscore<-function(mat,logr,logc)
{
score<-colSums((mat-rowMeans(mat)-matrix(colMeans(mat),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat))^2)/nrow(mat)
score
}

##Some helper funktions to calculate the CC score for node addition and inverse node addition

addrowscore<-function(mat,logr,logc)
{
score<-rowSums((mat-rowMeans(mat[,logc])-matrix(colMeans(mat[logr,]),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat[logr,logc]))^2)/ncol(mat[logr,logc])
score
}

iaddrowscore<-function(mat,logr,logc)
{
score<-rowSums((-mat+rowMeans(mat[,logc])-matrix(colMeans(mat[logr,]),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat[logr,logc]))^2)/ncol(mat[logr,logc])
score
}

addcolscore<-function(mat,logr,logc)
{
score<-colSums((mat-rowMeans(mat[,logc])-matrix(colMeans(mat[logr,]),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)+mean(mat[logr,logc]))^2)/nrow(mat[logr,logc])
score
}


# algorithm 1 from CC: Single Node Deletion

cc1<-function(mat,logr,logc,delta=1.5)
{
while(ccscore(mat[logr,logc])>delta)
{
di<-rowscore(mat[logr,logc])
dj<-colscore(mat[logr,logc])
mdi<-which.max(di)
mdj<-which.max(dj)

ifelse(di[mdi]>dj[mdj] ,logr[logr][mdi]<-FALSE ,logc[logc][mdj]<-FALSE)
if (!(sum(logr)>1 & sum(logc)>1))
break
}
ifelse(sum(logr)>1 & sum(logc)>1,ret<-list(logr,logc),ret<-list(0,warning(paste('No matirx with score smaller', delta,'found'))))
ret
}

# algorithm 2 from CC: Multiple Node Deletion
cc2<-function(mat,logr,logc,delta,alpha=1.5)
{
mdi<-1
mdj<-1
while((h<-ccscore(mat[logr,logc]))>delta & (sum(mdi)+sum(mdj))>0)
{

if(sum(logr)>100)
{
    di<-rowscore(mat[logr,logc])
    mdi<-di>(alpha*h)
    if(sum(mdi) < (sum(logr)-1))
    {
        logr[logr][mdi]<-FALSE
        h<-ccscore(mat[logr,logc])
    }
    else
    {
        print(warning(paste('Alpha', alpha,'to small!')))
        mdi <- 0
    }
}
else{mdi<-0}


if(sum(logc)>100)
{
    dj<-colscore(mat[logr,logc])
    mdj<-dj>(alpha*h)
    if(sum(mdj) < (sum(logc)-1))
    {
        logc[logc][mdj]<-FALSE
    }
    else
    {

        print(warning(paste('Alpha', alpha,'to small!')))
        mdi <- 0
    }
}
else{mdj<-0}
}

ret<-list(logr,logc)
ret
}

# algorithm 3 from CC:  Node Addition
cc3<-function(mat,logr,logc)
{
br<-1
ilogr<-rep(FALSE,length(logr))
while(br>0)
{
br1<-sum(logc)
br2<-sum(logr)
h<-ccscore(mat[logr,logc])
dj<-addcolscore(mat,logr,logc)

mdj<-dj<=h
logc[mdj]<-TRUE


h<-ccscore(mat[logr,logc])
di<-addrowscore(mat,logr,logc)
idi<-iaddrowscore(mat,logr,logc)

mdi<-di<=h
logr[mdi]<-TRUE
imdi<-idi<=h
mat[!(logr==imdi)&imdi]<- -mat[!(logr==imdi)&imdi]
logr[imdi]<-TRUE

br<-sum(logc)+sum(logr)-br1-br2
}
ret<-list(logr,logc)
ret

}



# Find biggest Bicluster:

bigcc<-function(mat,delta,alpha=1.5)
{
logr<-rep(TRUE,nrow(mat))
logc<-rep(TRUE,ncol(mat))
step1<-cc2(mat,logr,logc,delta,alpha)
step2<-cc1(mat,step1[[1]],step1[[2]],delta)
if(sum(step2[[1]])==0)
{ret<-list(0,warning(paste('Mo matrix with score smaller than', delta,'found')))
}
else{
ret<-cc3(mat,step2[[1]],step2[[2]])
}
ret
}


## Algorithm to find the number biggest bicluster.

ccbiclust<-function(mat,delta,alpha=1.5,number=100, rand=FALSE)
{
MYCALL <- match.call()
ma<-max(mat)
mi<-min(mat)
x<-matrix(FALSE,nrow=nrow(mat),ncol=number)
y<-matrix(FALSE,nrow=number,ncol=ncol(mat))
logr<-rep(TRUE,nrow(mat))
Stop <- FALSE
logr<-rep(TRUE,nrow(mat))
for(i in 1:number)
{
  if(sum(logr)<2)
  {
    Stop <- TRUE
    break
  }
  erg<-bigcc(mat[logr,],delta,alpha)
  if(sum(erg[[1]])==0)
  {
    Stop <- TRUE
    break
  }
  else
  {
    x[logr,i]<-erg[[1]]
    y[i,]<-erg[[2]]
    if(rand)
    {
        mat[erg[[1]],erg[[2]]]<-runif(sum(erg[[1]])*sum(erg[[2]]),mi,ma)
    }
    else
    {
        logr[logr][erg[[1]]] <- FALSE
    }
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


