cvmboost <- function(obj,ndat,nfold=10,figname=NULL,trace=FALSE){
       ### inspect coefficient path and AIC-based stopping criterion
       ### 10-fold cross-validation
       n <- ndat
       k <- nfold
       ntest <- floor(n / k)
      cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, n)), k - 1),
                         rep(0, n * k - (k - 1) * (n + ntest))), nrow = n)
       cvm <- cvrisk(obj, folds = cv10f)
       if(trace)
       print(cvm)
if(!is.null(figname)){
 postscript(figname, height = 6.9, width = 6.6,horizontal = FALSE, onefile = FALSE, print.it = FALSE)
 plot(cvm)
 dev.off()
}
        return(mstop(cvm))
}

