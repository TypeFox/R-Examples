"LjungBox" <-
function(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE){
     TestType <- "0"
    if (class(obj) == "ts" || class(obj) == "numeric" || class(obj) == 
        "matrix" || (class(obj)[1] == "mts" && class(obj)[2] == 
        "ts")) 
        TestType <- "1"
    if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == 
        "Arima" || class(obj) == "varest" || class(obj) == "FitAR" || 
        class(obj) == "FitFGN" || class(obj) == "garch" || 
        class(obj) == "fGARCH" || class(obj) == "list")
        TestType<-"2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, varest, FitAR, 
             FitFGN, garch, fGARCH, ts, numeric, matrix, (mts ts), or list")
     Maxlag <- max(lags)
     if (TestType=="1")
       res <- as.ts(obj)
     else{
             GetResid <- GetResiduals(obj)
             res <- GetResid$res
             order <- GetResid$order
       }
       if (SquaredQ){ 
         res <- res^2
         order <- 0
       }
       n <- NROW(res)
       k <- NCOL(res)
        df <- k^2*(lags-order)
        NegativeDF <- which(df<0)
        df[NegativeDF] <- 0
	     Accmat <- stats::acf(res, lag.max = Maxlag, plot = FALSE, type = "correlation")$acf
	     inveseR0 <- solve(Accmat[1,,])
             prodvec <- numeric(Maxlag)
	  for(l in 1:Maxlag){
            tvecR <- t(as.vector(Accmat[l+1,,]))
	        prodvec[l] <- 1/(n-l)*crossprod(t(tvecR),crossprod(t(kronecker(inveseR0,inveseR0)),t(tvecR)))
	  }
       Q <- n*(n+2)*cumsum(prodvec)
       STATISTIC <- Q[lags]
       PVAL <- 1 - stats::pchisq(STATISTIC,df)
       PVAL[NegativeDF] <- NA
       summary <- matrix(c(lags,STATISTIC,df,PVAL),ncol=4)
       dimnames(summary) <- list(rep("", length(STATISTIC)),c("Lags","Statistic","df","pvalue"))
    return(summary)
}
