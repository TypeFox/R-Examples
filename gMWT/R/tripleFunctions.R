# Version: 30-11-2012, Daniel Fischer

# This function calculates the test statistic for the triple type of our tests
# The outcommented lines are useful for the case, that we allow a general alternative
# formulated with a vector a. But here it only consumes calculation time...
getP.Rsub <- function(x=NULL,y=NULL,z=NULL,ST=NULL,EQ=NULL,ties=TRUE){
    if(is.null(ST)) ST <- createST(x,y,z)
    if(is.null(EQ)&&(ties==TRUE)) EQ <- createEQ(x,y,z)

    if(ties==TRUE){
      obsCounts <- sum(getSubMatrix(ST,1,2)%*%getSubMatrix(ST,2,3)) + 0.5*sum(getSubMatrix(ST,1,2)%*%getSubMatrix(EQ,2,3)) + 0.5*sum(getSubMatrix(EQ,1,2)%*%getSubMatrix(ST,2,3)) + (1/6)*sum(getSubMatrix(EQ,1,2)%*%getSubMatrix(EQ,2,3))
    } else {
      obsCounts <- sum(getSubMatrix(ST,1,2)%*%getSubMatrix(ST,2,3))
    }
    return(obsCounts/(ST$n1*ST$n2*ST$n3))
}

# # this function calculates the permutated values
perm.triple <- function(x,y,z,nper,algorithm){
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)

  result <- c()
  if(algorithm=="Rsubmat")
  {
#-----------
     ST <- createST(x,y,z)
     EQ <- createEQ(x,y,z)
     # Check if we have ties:
     ifelse(sum(EQ$mMatrix)>(Nx+Ny+Nz), ties <- TRUE, ties <- FALSE)
     STPerm <- ST
     EQPerm <- EQ
     
     if(ties==TRUE){
	for (i in 1:nper)
	{
	  perm <- sample(Nx+Ny+Nz)
	  STPerm$mMatrix <- ST$mMatrix[perm,perm]
	  EQPerm$mMatrix <- EQ$mMatrix[perm,perm]
 	
	  result[i] <- getP.Rsub(ST=STPerm,EQ=EQPerm,ties=TRUE)
	 }
      } else {
	for (i in 1:nper)
	{
	  perm <- sample(Nx+Ny+Nz)
	  STPerm$mMatrix <- ST$mMatrix[perm,perm]
	
	  result[i] <- getP.Rsub(ST=STPerm,ties=FALSE)
 	}
    }
  } else if (algorithm=="Cnaive"){
    for (i in 1:nper)
    {
      permvalues <- sample(c(x,y,z))
      result[i] <- getP.Cnaive(permvalues[1:Nx],permvalues[(Nx+1):(Nx+Ny)],permvalues[(Nx+Ny+1):(Nx+Ny+Nz)])
    }

  } else if (algorithm=="Rnaive"){
    for (i in 1:nper)
    {
      permvalues <- sample(c(x,y,z))
      result[i] <- getP.Rnaive(permvalues[1:Nx],permvalues[(Nx+1):(Nx+Ny)],permvalues[(Nx+Ny+1):(Nx+Ny+Nz)])
    }

  }else if(algorithm=="Csubmat"){
      result <- as.vector(getP.Csub(x,y,z,nper))
      #result <- result/(Nx*Ny*Nz)

  }else {
    stop("We can't do that!!!")
  }
  result
}