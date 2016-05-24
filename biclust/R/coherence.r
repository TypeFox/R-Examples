# Functions to compute different kinds of coherence measures on biclusters
# We consider measures to check sign coherent biclusters, 
# multiplicative and additive biclusters and constant biclusters, all of them
# by rows and/or columns.




# Equivalent function to diff() but that divides instead of substracting
# each row or column by the before/above row/column.
# x - matrix nxm to compute
# rowOrColumn - 1 makes divisions by rows, while 2 makes divisions by column.
#               Default is 1
# returns - a matrix (n-1)xm or nx(m-1) (depending on rowOrColumns) where row/
#           column i is the result of dividing row/column i+1 by row/column i
#           of input matrix.
# NOTE: divisions by zero are computed as Inf.

#Author: Rodrigo Santamaria
div=function(x,rowOrColumn=1)
  {
  n=dim(x)[1]
  m=dim(x)[2]

  if(rowOrColumn==1)
    {
    ret=matrix(NA,n,m)
    ret[1,]=rep(1,m)
    for(i in 2:n)
       {
       ret[i,]=x[i,]/x[(i-1),]
       }
    }
  else if(rowOrColumn==2)
    {
    ret=matrix(NA,n,m)
    ret[,1]=rep(1,n)
    for(j in 2:m)
       {
       ret[,j]=x[,j]/x[,(j-1)]
       }
    }
  ret
  }

# Adaptation of diff() to equal div() interface
# x - matrix nxm to compute
# rowOrColumn - 1 makes divisions by rows, while 2 makes divisions by column.
#               Default is 1
# returns - a matrix (n-1)xm or nx(m-1) (depending on rowOrColumns) where row/
#           column i is the result of dividing row/column i+1 by row/column i
#           of input matrix.
dif=function(x,rowOrColumn=1)
  {
  n=dim(x)[1]
  m=dim(x)[2]

  if(rowOrColumn==1)
    {
    ret=matrix(NA,n,m)
    ret[1,]=rep(0,m)
    for(i in 2:n)
       {
       ret[i,]=x[i,]-x[(i-1),]
       }
    }
  else if(rowOrColumn==2)
    {
    ret=matrix(NA,n,m)
    ret[,1]=rep(0,n)
    for(j in 2:m)
       {
       ret[,j]=x[,j]-x[,(j-1)]
       }
    }
  ret
  }

# As dif() and div() above, but only change of slope is extracted
# x - matrix nxm to compute
# rowOrColumn - 1 makes divisions by rows, while 2 makes divisions by column.
#               Default is 1
# returns - a matrix (n-1)xm or nx(m-1) (depending on rowOrColumns) where row/
#           column i is a vector on {-1,1} where aij is 1 if original row/column
#           (i-1) is smaller than row/column i, and -1 otherwise
sig=function(x,rowOrColumn=1)
  {
  n=dim(x)[1]
  m=dim(x)[2]

  if(rowOrColumn==1)
    {
    ret=matrix(NA,n,m)
    ret[1,]=rep(0,m)
    for(i in 2:n)
       {
       for(j in 1:m)
        {
        if(x[i,j]>x[(i-1),j]) ret[i,j]=1
        if(x[i,j]<x[(i-1),j]) ret[i,j]=-1
        if(x[i,j]==x[(i-1),j]) ret[i,j]=0
        }
       }
    }
  else if(rowOrColumn==2)
    {
    ret=matrix(NA,n,m)
    ret[,1]=rep(0,n)
    for(j in 2:m)
       {
       for(i in 1:n)
        {
        if(x[i,j]>x[i,(j-1)]) ret[i,j]=1
        if(x[i,j]<x[i,(j-1)]) ret[i,j]=-1
        if(x[i,j]==x[i,(j-1)]) ret[i,j]=0
        }
       }
    }
  ret
  }


#Computes variance of input matrix x by rows or columns as the
# sum of euclidean distances between all rows/columns, divided by 1/n(n-1),
# being n the number of rows or columns.
# Zero variance implies an homogeneous matrix by rows or columns
# BIC NOTE: This measure is only valid for biclusters based on constants values
# (by rows, columns, or both) or coherent evolutions based on grouping next values
# We hope it can be used by biclusters of coherent values (additive or multiplicative)
# by applying a previous differential treatment that not change variance measure
# in case of constant values.
# To do this, first apply differential treatmen with sig(), dif() or div()
# as follows
# Variance on genes -> variance(dif(M,2))
# Variance on conditions -> variance(dif(M),2)
variance=function(x,rowOrColumn=1)
  {
  n=dim(x)[1]
  m=dim(x)[2]
  ret=0
  
  if(rowOrColumn==1)         #rows
    {
    distan=as.matrix(dist(x))
    ret=(1/(n*(n-1)))*sum(dist(x))
    }
  else
    {
    if(rowOrColumn==2)  #cols
      {
      distan=as.matrix(dist(t(x)))
      ret=(1/(m*(m-1)))*sum(dist(t(x)))
      }
    else 
      {
      }
      
    }
    
  ret
  }

#Gives overall variance as mean of row and column variance() above
overallVariance=function(x)
  {
  rv=variance(x)
  cv=variance(x,2)
  n=dim(x)[1]
  m=dim(x)[2]
  
  ret=(n*rv+m*cv)/(n+m)
  ret
  }  

#Returns a dimxdim matrix P where pij=1/(n+1) where n is the number of 
#biclusters where rows/columns i and j appear grouped together. pij will be in
#[1/max,1] where 1 is achieved if they are never grouped and 1/max if they are
#grouped in all the biclusters of the result set.

  
constantVariance=function(x, resultSet, number, dimension="both")
  {
  rows=row(matrix(resultSet@RowxNumber[,number]))[resultSet@RowxNumber[,number]==T]
  cols=row(matrix(resultSet@NumberxCol[number,]))[resultSet@NumberxCol[number,]==T]
  A=x[rows,cols]

  if(dimension=="both")  
    {
    overallVariance(A)
    }
  else
    {
    if(dimension=="row")  variance(A,1)
    else  if(dimension=="col")  variance(A,2)
    }
  }
  
additiveVariance=function(x, resultSet, number, dimension="both")
  {
  if(dimension!="both" && dimension!="row" && dimension!="col")
    {
    stop("Argument dimension must be 'row', 'col' or 'both'")
    }
  rows=row(matrix(resultSet@RowxNumber[,number]))[resultSet@RowxNumber[,number]==T]
  cols=row(matrix(resultSet@NumberxCol[number,]))[resultSet@NumberxCol[number,]==T]
  A=x[rows,cols]

  if(dimension=="both")  
    {
    rv=variance(dif(A,2))
    cv=variance(dif(A),2)
    n=dim(A)[1]
    m=dim(A)[2]
  
    (n*rv+m*cv)/(n+m)
    }
  else
    {
    if(dimension=="row")  variance(dif(A,2))
    else  if(dimension=="col")  variance(dif(A),2)
    }
  }
  
multiplicativeVariance=function(x, resultSet, number, dimension="both")
  {
  if(dimension!="both" && dimension!="row" && dimension!="col")
    {
    stop("Argument dimension must be 'row', 'col' or 'both'")
    }
  rows=row(matrix(resultSet@RowxNumber[,number]))[resultSet@RowxNumber[,number]==T]
  cols=row(matrix(resultSet@NumberxCol[number,]))[resultSet@NumberxCol[number,]==T]
  A=x[rows,cols]

  if(dimension=="both")  
    {
    rv=variance(div(A,2))
    cv=variance(div(A),2)
    n=dim(A)[1]
    m=dim(A)[2]
  
    (n*rv+m*cv)/(n+m)
    }
  else
    {
    if(dimension=="row")  variance(div(A,2))
    else  if(dimension=="col")  variance(div(A),2)
    }
  }

signVariance=function(x, resultSet, number, dimension="both")
 {
  if(dimension!="both" && dimension!="row" && dimension!="col")
    {
    stop("Argument dimension must be 'row', 'col' or 'both'")
    }
  rows=row(matrix(resultSet@RowxNumber[,number]))[resultSet@RowxNumber[,number]==T]
  cols=row(matrix(resultSet@NumberxCol[number,]))[resultSet@NumberxCol[number,]==T]
  A=x[rows,cols]

  if(dimension=="both")  
    {
    rv=variance(sig(A,2))
    cv=variance(sig(A),2)
    n=dim(A)[1]
    m=dim(A)[2]
  
    (n*rv+m*cv)/(n+m)
    }
  else
    {
    if(dimension=="row")  variance(sig(A,2))
    else  if(dimension=="col")  variance(sig(A),2)
    }
  }
  
 
# Compute coherence measures for resultSet. resultSet must have a bicluster
# result syntax, concreetly $rows and $cols vectors for each bicluster, as
# returned by bimax() or plaid(). x is the matrix over which resultSet applies
# Returns row and column variance, additive, multiplicative and sign coherence
# for each biclusters, and overall quantities for all the result set:
#   resultSet$rv - row variance
#            $cv - col variance
#            $rac - row additive coherence
#            $cac - column additive coherence
#            $rmc - row multiplicative coherence
#            $cmc - column multiplicative coherence
#            $rsc - row sign coherence
#            $csc - column sign coherence
# 
coherenceAnalysis=function(x, resultSet)
  {
  overallRowVariance=0
  overallColVariance=0
  overallRowAddCoherence=0
  overallColAddCoherence=0
  overallRowMulCoherence=0
  overallColMulCoherence=0
  overallRowSignCoherence=0
  overallColSignCoherence=0
  N=0
  M=0
  ret=c()
  
  for(i in 1:length(resultSet))
    {
    rows=resultSet[[i]]$rows
    cols=resultSet[[i]]$cols
    A=x[rows,cols]
    ret[[i]]=list(rows=rows,cols=cols,
            rowVariance=variance(A,1), colVariance=variance(A,2), 
            rowAddCoherence=variance(dif(A,2)), colAddCoherence=variance(dif(A),2), 
            rowMulCoherence=variance(div(A,2)), colMulCoherence=variance(div(A),2),
            rowSignCoherence=variance(sig(A,2)), colSignCoherence=variance(sig(A),2) )
  
    N=N+length(rows)
    M=M+length(cols)
  overallRowVariance=overallRowVariance+length(rows)*ret[[i]]$rowVariance
  overallColVariance=overallColVariance+length(cols)*ret[[i]]$colVariance
  overallRowAddCoherence=overallRowAddCoherence+length(rows)*ret[[i]]$rowAddCoherence
  overallColAddCoherence=overallColAddCoherence+length(cols)*ret[[i]]$colAddCoherence
  overallRowMulCoherence=overallRowMulCoherence+length(rows)*ret[[i]]$rowMulCoherence
  overallColMulCoherence=overallColMulCoherence+length(cols)*ret[[i]]$colMulCoherence
  overallRowSignCoherence=overallRowSignCoherence+length(rows)*ret[[i]]$rowSignCoherence
  overallColSignCoherence=overallColSignCoherence+length(cols)*ret[[i]]$colSignCoherence
    }
    
  ret$rv=overallRowVariance/N
  ret$cv=overallColVariance/M
  ret$v=(ret$rv*N+ret$cv*M)/(N+M)
  ret$rac=overallRowAddCoherence/N
  ret$cac=overallColAddCoherence/M
  ret$ac=(ret$rac*N+ret$cac*M)/(N+M)
  ret$rmc=overallRowMulCoherence/N
  ret$cmc=overallColMulCoherence/M
  ret$mc=(ret$rmc*N+ret$cmc*M)/(N+M)
  ret$rsc=overallRowSignCoherence/N
  ret$csc=overallColSignCoherence/M
  ret$sc=(ret$rsc*N+ret$csc*M)/(N+M)
  ret
  }
  
#  bimaxCoherence=coherenceAnalysis(A,bicBimax)
#  plaidCoherence=coherenceAnalysis(A,bicPlaid)

