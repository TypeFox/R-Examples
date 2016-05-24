print.hiddenInvariantCausalPrediction <-
function(x,...){
    ns <- sum(x$maximinCoefficients!=0)
    wh <- which( x$maximinCoefficients!=0)
    cat(paste("\n Invariant Linear Causal Regression (with hidden variables) at level ", x$alpha,sep=""))

    if(ns>0) cat(paste("\n ",if(ns>1) "Variables: " else "Variable ",paste(x$colnames[wh[1:min(10,length(wh))]],collapse=", ")  ,if(length(wh)>10) paste("... (and ", length(wh)-10, "more) ",sep="") else ""," ", if(ns>1) "show" else "shows"," a significant causal effect",sep=""))
    sig <- apply(sign(x$ConfInt[2:1,,drop=FALSE]),2,function(x) prod(x))
    sigo <- sign(x$ConfInt[1,])
    dis <- rbind(x$ConfInt[2:1, ,drop=FALSE], sigo*apply(abs(x$ConfInt[2:1, ,drop=FALSE]),2,min) * (sig>=0),x$pvalues)
    colnames(dis) <- x$colnames
    rownames(dis) <- c( paste(" LOWER BOUND",sep=""),paste(" UPPER BOUND",sep=""),  " MAXIMIN EFFECT", " P-VALUE")
                                        #dis <- dis[ ,wh,drop=FALSE]
    cat("\n \n ")
    printCoefmat( t(dis), digits=3, signif.stars=TRUE,P.values=TRUE,has.Pvalue=TRUE,tst.ind=1:3,zap.ind=3,eps.Pvalue=10^(-9))
    
    cat("\n\n")
}
