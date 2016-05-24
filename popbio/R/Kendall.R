Kendall<-function(rates, grades=1000, maxvar=0.2,minvar=0.00001, maxmean=1, minmean=0.01)
{
   ## number of different years or time intervals and rate ids in dataset
   times<-length(unique(rates[,2]))
   classes<-unique(rates[,1])
   n<-length(classes)
   results <- matrix(rep(0,5*n), nrow=n,
               dimnames= list(classes, c("mean", "MLE.mean", "var", "MLE.var", "cor.MLE.var")))
   resultslik <-  matrix(rep(0,4*n), nrow=n,
               dimnames= list(classes, c("low.mean", "hi.mean", "low.var", "hi.var")))
   # makes all the sets of mean values to search over:
   means <- t(matrix(seq(minmean,maxmean,length = grades),grades,grades))
   # makes all the sets of variance values to search over:
   vars  <- matrix(seq(minvar,maxvar,length = grades),grades,grades)
   # using the means and variances to compute the a and b
   # (aa and bb) parameters of a beta distribution
   aa <- means * (((means * (1 - means)) / vars) - 1)
   bb <- (1 - means) * (((means * (1 - means)) / vars) - 1)

   # eliminate impossible combinations of mean and variance values
   aa[aa <= 0] <- NaN
   bb[bb <= 0] <- NaN
   
   for (clas in 1:n)     # loop through each rate or class
   {
      print(paste("Computing estimates for rate", classes[clas]), quote=FALSE)
      minrow   <- (clas - 1) * times + 1   # find the min and max rows of the data matrices to use
      maxrow   <-  clas * times
      dat      <- rates[minrow:maxrow,]    # fetch the data to use
      estmn    <- mean(dat[,4]/dat[,3],na.rm = TRUE)   # raw estimate of mean - excluding double zeros
      estvar   <- var(dat[,4]/dat[,3],na.rm = TRUE)    # raw est. of variance  - excluding double zeros
      # As is often done, use the negative of log liklihoods,
      # or -log likelihoods, in the following computations.
      loglikli <- matrix(0,grades,grades)              #initialize -loglikelihoods
      for (tt in 1:times)                              # calculate -loglikelihood for each year
      { 
         # this uses the -log-likelihood formula from Kendall:
         newlogL <-  -log(beta(aa + dat[tt,4],dat[tt,3]- dat[tt,4] + bb)/beta(aa,bb))  
         loglikli <- newlogL + loglikli                # add up the -log-likelihoods for each year
       
      } 
      # the lines below summarize the results
      minLL    <- min(loglikli,na.rm = TRUE)           # what is the best log-liklihood?
      minvars  <- apply(loglikli,2,min,na.rm = TRUE)   # finding best var
      ii       <- match(minvars,loglikli)
      minii    <- min(minvars,na.rm = TRUE)
      jj       <- which.min(minvars)
      MLEvar   <- vars[ii[jj]]
      minvars  <- apply(loglikli,1,min,na.rm = TRUE)   # finding best mean
      ii       <- match(minvars,loglikli)
      minii    <- min(minvars,na.rm = TRUE)
      jj       <- which.min(minvars)
      MLEmean  <- means[ii[jj]]
      clLL    <- loglikli - minLL                      # differences from best -LL
      hivar   <- max(vars[clLL<3.0],na.rm = TRUE)      # confidence limits for var
      lowvar  <- min(vars[clLL<3.0],na.rm = TRUE)
      himean  <- max(means[clLL<3.0],na.rm = TRUE)     # confidence limits for mean
      lowmean <- min(means[clLL<3.0],na.rm = TRUE)
      # perform correction for max. likelihood estimate of variance
      # (see Kendall 1998 for explanation)
      corrMLEvar <- MLEvar * times /(times - 1)
      corrlowvar <- lowvar * times /(times - 1)
      corrhivar  <- hivar * times /(times - 1)
      # storing the results
      results[clas,]  <-   c(estmn,MLEmean,estvar,MLEvar,corrMLEvar)
      resultslik[clas,] <- c(lowmean,himean,corrlowvar,corrhivar)
   }
   Kendall <- list(est = results,ci = resultslik)
   Kendall
}
