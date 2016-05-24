# Version: 30-11-2012, Daniel Fischer

# This function takes the given group vector and relabels it with 1 for the smallest
# value, 2 for the second smallest and so on
# Daniel, Tampere, 03-10-2012
# Function tested on 03-10-2012
  relabelGroups <- function(g){
    N <- length(g)
    gLables <- names(table(g))
    gReturn <- g 
    for(i in 1:length(gLables))
    {
      gReturn[g==gLables[i]] <- i
    }
    gReturn <- as.numeric(gReturn)
    return(gReturn)
  }
#----------------------------------------------------------------------------------------------------------------------------------------------

# This calculates the amount of possible probability estimatores, given a set of group labels, a type of PE and if the order matters:
# Daniel, Tampere, 03-10-2012
# Function tested on 03-10-2012
getNListItems <- function(Ngoi,type,order){
  res <- 0

  if(type=="pair")
  {
    if(order==TRUE)
    {
      res <- Ngoi*(Ngoi-1)/2
    } else {
      res <- Ngoi*(Ngoi-1)
    }
  } else if(type=="triple"){
    if(order==TRUE)
    {
      res <- Ngoi*(Ngoi-1)*(Ngoi-2)/6
    } else {
      res <- Ngoi*(Ngoi-1)*(Ngoi-2)
    }
  }

  return(res)
}
#----------------------------------------------------------------------------------------------------------------------------------------------


# The getComb function creates a matrix with possible testing groups for certain type and orders. options are the following:
# type="single"
#	Here we consider the group setting one group versus all other, each row is a different type of test in this setting,
#	and the first column contains the group label of the first testing group and the remaining ones the group labels of
#	the remaining groups.
# type="pairs"
#	All type of pair combinations of the given groups (result matrix contains 2 columns...). One additional column is possible,
#	the result can be ordered or not, in the order=F case both options 1 2 and 2 1 are given as rows, in the case order=T
#	we report only those cases where the first group is smaller than the second group.
# type="triple"
#	Here we give all triple combinations (hence the result matrix contains 3 columns). Again one additional argument is
#	possible, for the ordered case the first entry has to be smaller than the second one ond this one has to be smaller
#	than the third one. For order=F all combinations are provided,
# Daniel, Tampere 03-10-2012
# Function tested on 03-10-2012
getComb <- function(x,type,order){
  x <- sort(x)
  Nx <- length(x)
 
  if(type=="single"){
    result <- matrix(NA,ncol=Nx,nrow=Nx)
    for(i in 1:Nx)
    {
      result[i,1] <- x[i]
      result[i,2:Nx] <- x[-i]
    }
    result
  } else if(type=="pairs")
  {
    if(order==T)
    {
      result <- matrix(NA,ncol=2,nrow=Nx*(Nx-1)/2)
      row <- 1
      for(i in 1:(Nx-1))
      {
	for(j in (i+1):Nx)
	{
	  result[row,] <- c(x[i],x[j])
	  row <- row + 1
	}
      }
    } else {
      result <- matrix(NA,ncol=2,nrow=Nx*(Nx-1))
      row <- 1
      for(i in 1:(Nx-1))
      {
	for(j in (i+1):Nx)
	{
	  result[row,] <- c(x[i],x[j])
	  result[row + 1,] <- c(x[j],x[i])
	  row <- row + 2
	}
      }
    }
  } else if(type=="triple"){
    if(order==T)
    {
      result <- matrix(NA,ncol=3,nrow=Nx*(Nx-1)*(Nx-2)/6)
      row <- 1
      for(i in 1:(Nx-2))
      {
	for(j in (i+1):(Nx-1))
	{
	  for(k in (j+1):Nx)
	  {
	    result[row,] <- c(x[i],x[j],x[k])
	    row <- row + 1
	  }
	}
      }
    } else {
      result <- matrix(NA,ncol=3,nrow=Nx*(Nx-1)*(Nx-2))
      row <- 1
      for(i in 1:(Nx-2))
      {
	for(j in (i+1):(Nx-1))
	{
	  for(k in (j+1):Nx)
	  {
	    result[row,] <- c(x[i],x[j],x[k])
	    result[row+1,] <- c(x[i],x[k],x[j])
	    result[row+2,] <- c(x[j],x[i],x[k])
	    result[row+3,] <- c(x[j],x[k],x[i])
	    result[row+4,] <- c(x[k],x[i],x[j])
	    result[row+5,] <- c(x[k],x[j],x[i])
	    row <- row + 6
	  }
	}
      }
    }
  } else {
    result <- NULL
    stop("We do not have this type of combination in the offer today, sorry! Maybe some other day!?\n")
  }
 result
}
#----------------------------------------------------------------------------------------------------------------------------------------------

# This function calculates the U_{xy} hat statistic for two inputs and takes care of ties!
# UIT type test
getP.grid <- function(x,y){
 combinations <- expand.grid(x,y)
 # Cases for x < y
 tempSign <- (sign(combinations[,2] - combinations[,1]))
 result <- sum(tempSign==1)
 # Handle the ties
 result <- result + sum(tempSign==0)*0.5
 result/(length(x)*length(y))
}

getP.Rnaive <- function(x,y,z=NULL){
  result <- 0
  if(is.null(z))
  {
    for (i in 1:length(x))
    {
      for (j in 1:length(y))
      {
	if(x[i] < y[j]) result <- result + 1; 
	if(x[i] == y[j]) result <- result + 0.5;
      }
    }
    result <- result / (length(x)*length(y))
  } else {
    for (i in 1:length(x))
    {
      for (j in 1:length(y))
      {
	for(k in 1:length(z))
	{
	  if((x[i] < y[j]) & (y[j] < z[k])) result <- result + 1; 
	  if((x[i] == y[j]) & (y[j] < z[k])) result <- result + 0.5;
	  if((x[i] < y[j]) & (y[j] == z[k])) result <- result + 0.5;
	  if((x[i] == y[j]) & (y[j] == z[k])) result <- result + 1/6;
	}
      }
    }
    result <- result / (length(x)*length(y)*length(z))

  }
  result
}

# This function creates the indicator function matrix M. The design is as follows: It has #n1 + #n2 + #n3 columns and
# rows and binary entries. The observations are combined in the column and the rows as c(x,y,z), means the first n1 entries
# represents the observations from group z, and next n2 entries the ones from gtroup n2 and so on.
# The entries in the matrix are now the results of the pairwise comparison I(r[j]<c(i)), means, if the values in the row is
# smaller than the value in the column we put a 1, otherwise we put a 0.
# Additional the lengths of vectors x,y,z are in the output
createST <- function(x1,x2,x3){
  n1 <- length(x1)
  n2 <- length(x2)
  n3 <- length(x3)
  mMatrix.rows <- matrix(c(x1,x2,x3),ncol=(n1+n2+n3),nrow=n1+n2+n3,byrow=T)
  mMatrix.cols <- matrix(c(x1,x2,x3),ncol=(n1+n2+n3),nrow=n1+n2+n3,byrow=F)
  mMatrix <- mMatrix.rows-mMatrix.cols
  mMatrix <- (mMatrix>0)*1
  return(list(mMatrix=mMatrix,n1=n1,n2=n2,n3=n3))
}

createEQ <- function(x1,x2,x3){
  n1 <- length(x1)
  n2 <- length(x2)
  n3 <- length(x3)
  mMatrix.rows <- matrix(c(x1,x2,x3),ncol=(n1+n2+n3),nrow=n1+n2+n3,byrow=T)
  mMatrix.cols <- matrix(c(x1,x2,x3),ncol=(n1+n2+n3),nrow=n1+n2+n3,byrow=F)
  mMatrix <- mMatrix.rows-mMatrix.cols
  mMatrix <- (mMatrix==0)*1
  return(list(mMatrix=mMatrix,n1=n1,n2=n2,n3=n3))
}


# This function takes a mMatrix object (as created in the function above) and extractes the submatrix from it as specified.
# We have to give a mMatrix object and then two values i and j. This brings the submatrix for all indicator functions, that
# compare how often elements from group j (the later mentioned one) are greater than the ones mentioned in the first place.
# Hence, the submatrix has as many columns as elements are in the second mentioned group and as many rows as given in the first
# group
getSubMatrix <- function(mMatrix,i,j){
  if(i==1){
    rowLow <- 1
    rowUp <- mMatrix$n1
  } else if(i==2){
    rowLow <- mMatrix$n1+1
    rowUp <- mMatrix$n1+mMatrix$n2
  } else if (i==3){
    rowLow <- mMatrix$n1+mMatrix$n2 +1
    rowUp <- mMatrix$n1+mMatrix$n2+ mMatrix$n3
  }
  if(j==1){
    colLow <- 1
    colUp <- mMatrix$n1
  } else if(j==2){
    colLow <- mMatrix$n1+1
    colUp <- mMatrix$n1+mMatrix$n2
  } else if (j==3){
    colLow <- mMatrix$n1+mMatrix$n2 +1
    colUp <- mMatrix$n1+mMatrix$n2+ mMatrix$n3
  }
  return(mMatrix$mMatrix[rowLow:rowUp,colLow:colUp])
}
#----------------------------------------------------------------------------------------------------------------------------

# This function takes an estProb element and a value and returns the probability vector
getSingleProb <- function(x,i){
  if(is.matrix(x$probs))
  {
    tempPos <- rownames(x$probs) == paste("P(",i,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos,]))
  } else {
    tempPos <- names(x$probs) == paste("P(",i,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos]))
  }
 return(result)
}

getPairProb <- function(x,i,j){
  if(is.matrix(x$probs))
  {
    tempPos <- rownames(x$probs) == paste("P(",i,"<",j,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos,]))
  } else {
    tempPos <- names(x$probs) == paste("P(",i,"<",j,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos]))
  }
 return(result)
}

getTripleProb <- function(x,i,j,k){
  if(is.matrix(x$probs))
  {
    tempPos <- rownames(x$probs) == paste("P(",i,"<",j,"<",k,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos,]))
  } else {
    tempPos <- names(x$probs) == paste("P(",i,"<",j,"<",k,")",sep="")
    result <- as.vector(unlist(x$probs[tempPos]))
  }
 return(result)

}

# Combinations for the plotting of pair probability estimators
getPairComb <- function(x,order){
  Nx <- length(x)
  combs <- Nx*(Nx-1)/2
 # totalCombs <- combs*(combs-1)/2
    
  pairs <- getComb(x,type="pairs",order=order)
  rowCombs <- getComb(1:nrow(pairs),type="pairs",order=TRUE)

  result <- matrix(NA,ncol=4,nrow=nrow(rowCombs))

  for(i in 1:nrow(rowCombs))
  {
    result[i,] <- c(pairs[rowCombs[i,1],]  , pairs[rowCombs[i,2],])
  }
 result
}

# Combinations for the plotting of triple probability estimators
# Combinations for the plotting of triple probability estimators
getTripleComb <- function(x,order){
   Nx <- length(x)
   combs <- Nx*(Nx-1)*(Nx-2)
       
   triples <- getComb(x,type="triple",order=order)
   rowCombs <- getComb(1:nrow(triples),type="pairs",order=TRUE)

   result <- matrix(NA,ncol=6,nrow=nrow(rowCombs))

   for(i in 1:nrow(rowCombs))
   {
     result[i,] <- c(triples[rowCombs[i,1],]  , triples[rowCombs[i,2],])
   }
  result
}

# This function is needed for the P_t case. The two other cases are done in C, this one here only combines the pairs to get  the
# Probability estimator for P_t
oneVsOther <- function(X,g,t,goi){
   Nt <- sum(g==t)
   N <- length(g)
   result <- 0
   tprime <- unique(g[g!=t])
   for(i in 1:length(tprime))
   {
      result <- result + sum(g==tprime[i])*getP.Cnaive(X[g==t],X[g==tprime[i]])
   }
   result/(N-Nt)
}