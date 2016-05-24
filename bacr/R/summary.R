summary.bacr <- function(object, ...){
    output = object$output
    summ = list(posterior.mean=mean(object$ACE), CI=quantile(object$ACE, c(0.025, 0.975)))
    vv = colMeans(object$MY)
    names(vv) = object$predictorsY
    vv = vv[1:(length(vv)-1)]
    summ$PIP = vv
    
    cat("\nBAC objects:\n\n")
    ww = c(paste(signif(summ$posterior.mean,2)), paste("(",signif(summ$CI[1],2), ", ", signif(summ$CI[2],2),")",sep=""))
    names(ww) = c("posterior mean", "95% posterior interval")
    
    cat("Exposure effect estimate:\n")
    print(ww, print.gap = 5, quote = FALSE)
    
    qq = sort(summ$PIP, decreasing = T)
    qq = as.matrix(qq[qq>0.5])
    if (nrow(qq)>0){
        cat("\n\nCovariates with posterior inclusion probability > 0.5:\n")
        colnames(qq) = "posterior inclusion probability"
        print(qq)
    } else{
        cat("\n\nThere is no covariate with posterior inclusion probability > 0.5:\n")
    }
    
    invisible(summ)
}