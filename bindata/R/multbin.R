### convert a real valued array to a binary array by thresholding at zero

ra2ba <- function(x) {

  retval <- as.numeric(x>0)
  dim(retval) <- dim(x)
  retval
}


### generate multivariate binary data by direct conversion
### from multivariate normals   

rmvbin <- function (n, margprob,
                      commonprob=diag(margprob),
                      bincorr=diag(length(margprob)),
                      sigma=diag(length(margprob)),                      
                      colnames=NULL, simulvals=NULL) {

  if(missing(sigma))
    {
      if(!missing(commonprob))
        {
          if (missing(margprob))
            margprob <- diag(commonprob)
          sigma <- commonprob2sigma(commonprob, simulvals)
        }
      else if(!missing(bincorr))
        {
          commonprob <- bincorr2commonprob(margprob, bincorr)
          sigma <- commonprob2sigma(commonprob, simulvals)
        }
    }
  else if (any(eigen(sigma)$values<0))
    stop ("Sigma is not positive definite.")
  
  retval <- rmvnorm(n, qnorm(margprob), as.matrix(sigma))
  retval <- ra2ba(retval)
  dimnames(retval) <- list(NULL, colnames)
  retval
}


### convert binary correlation matrix to matrix of joint probabilities

bincorr2commonprob <- function(margprob, bincorr) {

  retval <- 0 * bincorr

  for(k in 1:ncol(retval)){
    for(l in 1:ncol(retval)){
      retval[k,l] <- bincorr[k,l] *
        sqrt(margprob[k]*(1-margprob[k])*margprob[l]*(1-margprob[l])) +
        margprob[k]*margprob[l]
    }
  }

  retval
}
      
### convert matrix of joint probabilities to covariance matrix of normal
### distribution 

commonprob2sigma <- function(commonprob, simulvals=NULL) {

  if(is.null(simulvals)){
      simulvals <- SimulVals
  }
      
  
  margprob <- diag(commonprob)
  
  if(!(check <- check.commonprob(commonprob))){
    cat(attr(check, "message"), sep="\n")
    stop("Matrix commonprob not admissible.")
  }
  
  sigma <- diag(nrow(commonprob))
  
  for(m in 1:(ncol(commonprob)-1)){
    for(n in (m+1):nrow(commonprob)){
      x <- cbind(margprob[m], margprob[n],
                 as.numeric(dimnames(simulvals)[[3]]))
      y <- interpolate(x, simulvals)
      f <- approxfun(y, x[,3])
      sigma[m,n] <- sigma[n,m] <- f(commonprob[m,n])
    }
  }
  if(any(is.na(sigma)))
    stop("Extrapolation occurred ... margprob and commonprob not compatible?")
  if (any(eigen(sigma)$values < 0))
    {
      cat ("Warning: Resulting covariance matrix is not positive definite.\n")
      cat ("         Smallest eigenvalue equals", min(eigen(sigma)$values),
           ".\n")
      cat ("         Please check whether the results are still useful.\n")
    }
  
  sigma
}




### compute the conditional probabilities P(x_i|x_j)

condprob <- function (x) {

  x <- as.matrix(x)
  nc <- ncol(x)

  retval <- matrix(0, nrow=nc, ncol=nc)
  for(k in 1:nc){
    retval [k,] <- apply(x[x[,k]!=0,], 2, mean)
  }
  dimnames(retval) <- list(colnames(x), colnames(x))
  retval
}
      

### check the matrix of common probabilities

check.commonprob <- function (commonprob)
{
  retval <- TRUE
  message <- character(0)
  nm <- 0
  
  if ((any(commonprob < 0)) || (any(commonprob > 1))){
    retval <- FALSE
    message[nm<-nm+1] <- "Not all probabilities are between 0 and 1."
  }
  
  n <- dim(commonprob)[1]
  if (n != dim(commonprob)[2]){
    retval <- FALSE
    message[nm<-nm+1] <- "Matrix of common probabilities is not quadratic."
  }
  
  ## check pairwise conditions
  for (i in 1:(n-1)){
    for (j in 2:n)
      {
        ul <- min(commonprob[i,i], commonprob[j,j])
        ll <- max(commonprob[i,i]+commonprob[j,j]-1, 0)
        if ((commonprob[i,j] > ul) || (commonprob[i,j] < ll))
          {
            retval <- FALSE
            message[nm<-nm+1] <- 
              paste("Error in Element (",i,",",j,
                    "): Admissible values are in [",
                    ll,",",ul,"].")
          }
      }
  }
  
  ## check triple conditions
  if (n > 2)
    for (i in 1:(n-2))
      for (j in (i+1):(n-1))
        for (k in (j+1):n)
          {
            l <- commonprob[i,i]+commonprob[j,j]+commonprob[k,k]-1
            if (commonprob[i,j]+commonprob[i,k]+commonprob[j,k] < l)
              {
                retval <- FALSE
                message[nm<-nm+1] <-
                  paste("The sum of the common probabilities of",i,",",
                        j,",",k,"must be at least",l,".")
              }
          }

  attr(retval, "message") <- message
  retval
}

