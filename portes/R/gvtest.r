"gvtest" <-
function(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE,Kernel=FALSE){
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
     k <- NCOL(res)
     n <- NROW(res)
    Det <- numeric(length(lags))
    if (Kernel==FALSE){
    mat <- ToeplitzBlock(res, lag.max=max(lags))
    }
    else {

      mat <- ToeplitzBlock(res, lag.max=max(lags),Kernel=TRUE)
    }

    for (i in 1:length(lags))
    Det[i] <- (-3*n/(2*lags[i]+1))*log(det(mat[(1:((lags[i] +1 ) * k)), (1:((lags[i] + 1) * k))]))
    df <- k^2*(1.5*lags*(lags+1)/(2*lags+1)-order)
    NegativeDF <- which(df<0)
    df[NegativeDF] <- 0
    PVAL <- 1 - stats::pchisq(Det,df)
    PVAL[NegativeDF] <- NA
    summary <- matrix(c(lags,Det,df,PVAL),ncol=4)
    dimnames(summary) <- list(rep("", length(Det)),c("Lags","Statistic","df","pvalue"))
  return(summary)
}
