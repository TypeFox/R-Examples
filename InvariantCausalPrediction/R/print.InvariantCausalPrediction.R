print.InvariantCausalPrediction <-
function(x,...){
    ns <- sum(x$maximinCoefficients!=0)
    wh <- which( x$maximinCoefficients!=0)
    cat(paste("\n Invariant Linear Causal", if(x$factor)" Classification " else " Regression ", "at level ", x$alpha, " (including multiplicity correction for the number of variables)",sep=""))
    if(x$modelReject){
        cat(paste("\n Model has been rejected at the chosen level ",x$gof,", that is no subset of variables leads to invariance across the environments. This can be for example due to  presence of \n (a) non-linearities or \n (b) hidden variables or \n (c) interventions on the target variable. \n",sep=""))
        if(x$stopIfEmpty){
            cat("\n In this run, option 'stopIfEmpty' was set to TRUE so not all sets of variables have maybe been tested; rerun with option set to FALSE to get definite answer whether model is rejected")
        }
          cat(paste("\n We will try to extend the functionality soon to allow non-linear models and address issue (a) [non-linearity], which currently leads to rejection of the linear model."))
          cat(paste("\n If the reason might be related to issue (b) [presence of hidden variables], one can use function hiddenICP which allows for hidden variables."))
          if(x$noEnv>2) cat(paste("\n If the reason might be related to issue (c) [interventions on the target in some environments], could repeat analysis while using only two of the ", x$noEnv," environments in one analysis and then taking the union of causal effects over all pairs of environments (interventions on the target in a given pair of environments will in general yield an empty set of estimated causal variables) ",sep=""))
    }else{
        if(ns>0) cat(paste("\n ",if(ns>1) "Variables: " else "Variable ",paste(x$colnames[wh[1:min(10,length(wh))]],collapse=", ")  ,if(length(wh)>10) paste("... (and ", length(wh)-10, "more) ",sep="") else ""," ", if(ns>1) "show" else "shows"," a significant causal effect",sep=""))
        sig <- apply(sign(x$ConfInt[2:1, ,drop=FALSE]),2,function(x) prod(x))
        sigo <- sign(x$ConfInt[1,])
        dis <- rbind(x$ConfInt[2:1, ,drop=FALSE], sigo*apply(abs(x$ConfInt[2:1, ,drop=FALSE]),2,min) * (sig>=0),x$pvalues)
        colnames(dis) <- x$colnames
        rownames(dis) <- c( paste(" LOWER BOUND",sep=""),paste(" UPPER BOUND",sep=""),  " MAXIMIN EFFECT", " P-VALUE")
                                        #dis <- dis[ ,wh,drop=FALSE]
        cat("\n \n ")
        printCoefmat( t(dis), digits=3, signif.stars=TRUE,P.values=TRUE,has.Pvalue=TRUE,tst.ind=1:3,zap.ind=3,eps.Pvalue=10^(-9))
        
    }
    cat("\n\n")
}
