#################################################################
#
# stutil.R
#
#################
# stepp utility #
#################
library(methods)
library(car)
#library(cmprsk)

setGeneric("estimate", function(.Object,...)
	standardGeneric("estimate"))

setGeneric("test", function(.Object,...)
	standardGeneric("test"))	

#
# internal routine to generate covariance matrix, pvalue
#
ssigma <- function(imatrix)
{
      sigma <- var(imatrix)

#     catch the error if we encounter a singularity problem when inverting the matrix sigma
#
      tryCatch({
      	sigmainv <- solve(sigma)
	      }, error = function(ex){
		cat(geterrmessage())
		cat("\n")
		set.seed<-4593432
		stop("STEPP has encountered an unexpected problem. Please retry. If problem persists, try different values of minpatspop(r1) and patspop(r2).")
	  	}
        )
     
	return (list(sigma=sigma, sigmainv=sigmainv))
}

#
# internal routine to generate the permutation pvalue
# compute the chisq or homogeneous association statistics and generate the pvalue
#
ppv <- function(imatrix, sigmainv, estarray, est, noPerms)
{
      perm <- rep(NA,noPerms)
      for (i in 1:noPerms){
        temp <- matrix(imatrix[i,],ncol=1)
        perm[i] <- t(temp)%*%sigmainv%*%temp
      }
      obs <- matrix(estarray - est, ncol=1)
      obs1 <- c(t(obs)%*%sigmainv%*%obs)
      pvalue <- sum(ifelse(perm > obs1, 1, 0))/noPerms
	signif(pvalue, 6)

      return(pvalue)
}

# compute the supremum statistics and generate the pvalue
#
ppv2 <- function(imatrix, estarray, est, noPerms, debug=FALSE)
{
#
      sigma <- sqrt(diag(var(imatrix)))
      stdDifferences <- t(apply(imatrix, 1, function(x) x/sigma))
      tPerm <- apply(abs(stdDifferences), 1, max)
 
      obsDifferences <- t(matrix(estarray - est, ncol=1))
      stdObsDifferences <- apply(obsDifferences, 1, function(x) x/sigma)
      tObs <- apply(abs(stdObsDifferences), 2, max)
      pvalue <- sum(ifelse(tPerm > tObs, 1, 0))/noPerms
      signif(pvalue, 6)

	return(pvalue)
}

