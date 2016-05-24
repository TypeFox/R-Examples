#' @S3method summary Roc
summary.Roc <- function(object,digits=2,print.it=TRUE,...){
    if (print.it) cat("Receiver operating characteristic\n")
    if (!is.null(object$method) && print.it) print(object$method)
    if (print.it){
        tabR <- table(object$Response)
        textR <- sapply(1:2,function(i){
            paste("'",names(tabR)[i],"' (n=",tabR[i],")",sep="")
        })
        cat("\nResponse: ",paste(textR,collapse="\t"),"\n\n")
    }
    if (length(object$models)>0){
        # ------------------------Area under the curve------------------------
        modelsAuc <- do.call("rbind",lapply(object$Auc,function(m){round(100*unlist(m),digits=digits)}))
        rownames(modelsAuc) <- names(object$models)
        xxx <- names(object$Auc[[1]])
        names(xxx) <- xxx
        xxx[xxx=="Auc"] <- object$method$name
        colnames(modelsAuc) <- xxx
        if (!is.null(object$BootcvRocList)){
            modelsAuc.ci <- lapply(object$BootcvRocList, function(y)quantile(sapply(y,function(z)Auc.default(z$Sensitivity,z$Specificity)), c(0.025,0.975)))
            ci <- round(100*do.call("rbind",modelsAuc.ci),digits)
            colnames(ci) <- paste("bootcv.",colnames(ci),sep="")
            modelsAuc <- cbind(modelsAuc,ci)
        }

        # ----------------------------Brier score----------------------------
        modelsBS <- do.call("rbind",lapply(object$Brier,function(m){round(100*unlist(m),digits=digits)}))
        rownames(modelsBS) <- names(object$models)
        xxx <- names(object$Brier[[1]])
        names(xxx) <- xxx
        xxx[xxx=="BS"] <- object$method$name
        colnames(modelsBS) <- xxx
        if (print.it) {
            if (object$method$name=="noPlan"){
                cat("Area under the ROC curve (AUC, higher better)\nBrier score (Brier, lower better)\n\n")
                if (NROW(modelsAuc)==1)
                    tab <- cbind(modelsAuc[order(-modelsAuc[,1,drop=TRUE]),,drop=FALSE],
                                 modelsBS[order(modelsBS[,1,drop=TRUE]),,drop=FALSE])
                else
                    tab <- cbind(modelsAuc[order(-modelsAuc[,1,drop=TRUE]),,drop=FALSE],
                                 modelsBS[order(modelsBS[,1,drop=TRUE]),,drop=FALSE])
                colnames(tab) <- c("AUC","Brier")
                print(tab)
                cat("\nEither newdata or apparent (learn data) performance.\n")
            } else{
                cat("Area under the ROC curve (AUC, higher better):\n")
                tab.auc <- apply(as.data.frame(modelsAuc[order(-modelsAuc[,1,drop=TRUE]),,drop=FALSE]),2,round,digits=digits)
                if (is.null(dim(tab.auc)))
                    names(tab.auc) <- gsub("AppAuc","apparent",names(tab.auc))
                else
                    colnames(tab.auc) <- gsub("AppAuc","apparent",colnames(tab.auc))
                print(tab.auc,quote=FALSE)
                cat("\nBrier score (Brier, lower better):\n")
                tab.brier <- apply(as.data.frame(modelsBS[order(modelsBS[,1,drop=TRUE]),,drop=FALSE]),2,round,digits=digits)
                if (is.null(dim(tab.auc)))
                    names(tab.brier) <- gsub("AppBS","apparent",names(tab.brier))
                else
                    colnames(tab.brier) <- gsub("AppBS","apparent",colnames(tab.brier))
                print(tab.brier, quote=FALSE)
            }
        }
        out <- list(Auc=modelsAuc,Brier=modelsBS)
        class(out) <- "summary.Roc"
        invisible(out)
    }
}

#' @S3method print summary.Roc
print.summary.Roc <- function(x,digits=2,...){
  lapply(x,function(x){
    cat("\n")
    print(apply(as.data.frame(x),2,round,digits=digits),quote=FALSE)
  })
}

##   cat("\n\n")
##   if (!is.null(object$PPV)){
##     PV(x)


## PV <- function(object)
## {
##   if (!is.null(object$confint) && object$confint==TRUE){
##     cat("Confidence method: ",object$confint.method)
##     cat("\n\n")
##     prSens <- object$Roc$Sensitivity[who]
##     prSpec <- object$Roc$Specificity[who]
##     prPPV <- object$Roc$PPV[who]
##     prNPV <- object$Roc$NPV[who]
##     myDigits <- if (percent==TRUE) digits-2 else digits
##     myFactor <- if (percent==TRUE) 100 else 1 
##     out <- do.call("cbind",list(Sens=round(myFactor*prSens,myDigits),"(CI.Sens)"=MgFormCi(lower=object$CI.Sens[who,"Lower"],upper=object$CI.Sens[who,"Upper"]),Spec=round(myFactor*prSpec,myDigits),"(CI.Spec)"=MgFormCi(lower=object$CI.Spec[who,"Lower"],upper=object$CI.Spec[who,"Upper"]),PPV=round(myFactor*prPPV,myDigits),"(CI.PPV)"=MgFormCi(lower=object$CI.PPV[who,"Lower"],upper=object$CI.PPV[who,"Upper"]),NPV=round(myFactor*prNPV,myDigits),"(CI.NPV)"=MgFormCi(lower=object$CI.NPV[who,"Lower"],upper=object$CI.NPV[who,"Upper"])))
##     LRplus <- prSens/(1-prSpec)
##     LRminus <- (1-prSens)/prSpec
##     if (length(breaks)==NROW(out)-1)
##       rownames(out) <- c(round(breaks,digits),"--")
##     else
##       rownames(out) <- round(breaks,digits)
##     cat("\n\nSens=Sensitivity\nSpec=Specificity\nPPV=positive predictive value\nNPV=negative predictive value\n\n")
##     print(out,quote=FALSE)
##     cat("\n\nLR+=positive likelihood ratio\nLR-=negative likelihood ratio\n\n")
##     print(cbind(rownames(out),"LR+"=round(LRplus,digits),"LR-"=round(LRminus,digits)),quote=FALSE)
##   }
##   else{
##     out <- do.call("cbind",list(Sens=object$Roc$Sensitivity[who],Spec=object$Roc$Specificity[who],PPV=object$Roc$PPV[who],NPV=object$Roc$NPV[who]))
##     rownames(out) <- round(breaks,digits)
##     cat("\n\nSens=Sensitivity\nSpec=Specificity\nPPV=positive predictive value\nNPV=negative predictive value\n\n")
##     print(100*out,digits=digits)
##   }
## }
