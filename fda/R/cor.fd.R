cor.fd <- function(evalarg1, fdobj1, evalarg2=evalarg1, fdobj2=fdobj1)
{
  #  compute correlation matrix / matrices for functional observations

  #  Last modified 25 February 2007 by Spencer Graves
##
## 1.  Compute var1 = bivariate data object for variance(fdobj1)
##
  var1 <- var.fd(fdobj1)
##
## 2.  Evaluate var1 at evalarg1 
##
  evalVar1 <- eval.bifd(evalarg1, evalarg1, var1)
##
## 3.  If missing(fdobj2) convert evalVar1 to correlations 
##
  {
    if(missing(fdobj2)){
      dimV1 <- dim(evalVar1)
      ndV1 <- length(dimV1)
      {
        if(ndV1<3){
          s1 <- sqrt(diag(evalVar1))
          return(evalVar1/outer(s1, s1))
        }
        else{
          if(dimV1[3] != 1)
            stop("Bug in cor.fd:  Programmed only for ",
                 "matrices or 4-d arrays with dim[3]==1.  oops.")
          dim1 <- dim(fdobj1$coefs)
          nVars <- dim1[3]
#         The following identifies the levels of evalVar1
#         containing the covariance matrix of each variable with itself           
          evalV.diag <- choose(2:(nVars+1), 2)
#         Compute the standard deviation vector for each variable 
          nPts <- length(evalarg1)
          s1 <- array(NA, dim=c(nPts, nVars))
          for(i in 1:nVars)
            s1[, i] <- sqrt(diag(evalVar1[,,1,evalV.diag[i]]))
#         Now compute the correlations 
          cor1 <- evalVar1
          m <- 0
          for(i in 1:nVars)for(j in 1:i){
            m <- m+1
            cor1[,,1,m] <- (evalVar1[,,1,m] / outer(s1[, i], s1[, j]))
          }
#         end for i in 1:nVars and j in 1:i
          return(cor1)
        }
#       end if...else nV1>2
      }
    }
    else {
##
## 4.  fdobj2 was also provided 
##
#  4.1  var.df(fdobj2)       
      var2 <- var.fd(fdobj2)
#  4.2.  Evalate at evalarg2
      evalVar2 <- eval.bifd(evalarg2, evalarg2, var2)
#  4.3.  var12 cross covariance
#*** If fdobj1 or fdobj2 are multivariate, var.fd will complain.        
      var12 <- var.fd(fdobj1, fdobj2)
#  4.4.  Evaluate the cross covariances           
      evalVar12 <- eval.bifd(evalarg1, evalarg2, var12)
#  4.5.  Convert evalVar12 to correlations 
      s1 <- sqrt(diag(evalVar1))
      s2 <- sqrt(diag(evalVar2))
      return(evalVar12 / outer(s1, s2))
    }
  }
##
## Done:  Nothing should get this far;
## all should have returned earlier,
## either via univariate or multivariate
## with fdobj1 only 
## or with both fdobj1 and fdobj2
##   
}
