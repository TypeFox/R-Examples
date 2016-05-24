# wrapper for CKricoefs...
OKrig <- function (x, Y)
{
   xy <- x
   xy$y <- Y
   nuniquerows <- nrow(unique(x))
   out <- CKrigcoefs(xy, nuniquerows, covfnparamA=blackbox.getOption("CovFnParam"),
                      lambdaA=blackbox.getOption("lambdaEst")) ## newCSmooth + Krig_coef_Wrapper + deleteCSmooth
   out$nuniquerows <- nuniquerows
   dim(out$c) <-c(length(out$c),1L)
   out$d <- out$d[1L] ## $d is full length but only first elemnt is non-zero for Ordinary Kriging
   dim(out$d) <-c(1L,1L)
   out <- c(out,list(x=x,y=Y))
   class(out) <- c("OKrig",class(out)) ## FR->FR not sur it is still used
   return(out)
}

