spectral=function(x,normalization="log", numberOfEigenvalues=3, 
            minr=2, minc=2, withinVar=1)
  {
  A=x
  n=dim(A)[1]
  m=dim(A)[2]
  discardFirstEigenvalue=T
  
  #1) Normalization------
  if(normalization=="log")
    {
    if(min(A)<1)  A= A+1+(abs(min(A)))
    K=logt(A)
    discardFirstEigenvalue=F
    }
  else
    {
    if(numberOfEigenvalues==1) 
      {
      warning("Only 1 eigenvalue selected with a normalization that gives first eigenvalue as background. Increasing to 2 eigenvalues")
      numberOfEigenvalues=2;
      }
      
    Ab=A+abs(min(A))
    if(normalization=="irrc")                   K=irrc(Ab)
    else if(normalization=="bistochastization") K=bistochastization(Ab)
    }
    
  #2) SVD Decomposition
  desc=svd(K)

  #3) Vector processing: partitive clustering and reordering
#  result=postprocess(desc,maxeigen=numberOfEigenvalues,minCG=max(2,floor(n/400)),
#          maxCG=min(n/minr,100),minCE=max(2,floor(m/20)),maxCE=m/minc)
   result=postprocess(desc,maxeigen=numberOfEigenvalues,minCG=2,
          maxCG=min(n/minr,100),minCE=2,maxCE=m/minc)

  #4) Taking biclusters of all possible eigenvector combinations
  srows=list()
  scols=list()
  init=1
  if(discardFirstEigenvalue) init=2
  for(k in init:numberOfEigenvalues) 
    {
    init2=k
    for(l in init2:numberOfEigenvalues)
      {
      numGenes=result$numgene[k]
      numExpr=result$numexpr[l]

      for(i in 1:numGenes)
        {
        for(j in 1:numExpr)
          {
          srows=c(srows,list(row(A)[result$eigengenecluster[k,]==i,1]))
          scols=c(scols,list(col(A)[1,result$eigenexprcluster[l,]==j]))
          }
        }
      }
    }

  #5) Discarding non-relevant biclusters
  srowsOK=list()
  scolsOK=list()
  ss=c()
  ret=list()
  for(i in 1:length(srows))
    {
    ncols=length(scols[[i]])
    nrows=length(srows[[i]])
    wv=withinVar(A[srows[[i]],scols[[i]]],ncols,nrows)
    if(wv<withinVar && ncols>minc && nrows>minr)
      {
       srowsOK=c(srowsOK, list(srows[[i]]))
       scolsOK=c(scolsOK, list(scols[[i]]))
       ss=c(ss,wv)
      }
    }

    if(length(srowsOK)==0)
      {
      warning("No biclusters found")
      return(BiclustResult(as.list(match.call()),matrix(NA,1,1),matrix(NA,1,1),0,list(0)))
      }
    else
      {
      rowxnumber=matrix(F, nrow=dim(A)[1], ncol=length(srowsOK))
      colxnumber=matrix(F, nrow=dim(A)[2], ncol=length(srowsOK))
    
      for(i in 1:length(srowsOK))
        {
        temp=rep(FALSE,dim(A)[1])
        temp[srowsOK[[i]]]=T
        rowxnumber[,i]=temp
       
        temp=rep(FALSE,dim(A)[2])
        temp[scolsOK[[i]]]=T
        colxnumber[,i]=temp
        }
      return(BiclustResult(as.list(match.call()),rowxnumber,t(colxnumber),length(srowsOK),list(0)))
      }
  }