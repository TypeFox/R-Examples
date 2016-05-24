extractcoef <- function(x){
    if (x$convergence == 1) {
        ms <- x$coefficients
        mnames <- rownames(ms)
        nrows <- dim(ms)[1]
        newnames <- paste0(mnames,c(rep("",nrows),rep("_se",nrows),
                                    rep("_t",nrows),rep("_pval",nrows)))
        aux_res <- as.vector(ms)
        names(aux_res) <- newnames
        res <- data.frame(t(aux_res))
        res$obs <- nrow(x$data)
        res$rsquare <- x$rsquare
        res$modelfit <- x$modelfit
        res$aic <- x$aic
        res$convergence <- x$convergence
        res$plateid <- x$data$plateid[1]
        res$fct <- x$fct
        res$bkg_mean <- x$bkg_mean
        res$bkg_method <- x$bkg_method
    } else {
        res <- data.frame(obs=nrow(x$data), convergence=x$convergence, 
                        plateid=x$data$plateid[1], fct=x$fct, 
                        bkg_mean=x$bkg_mean, bkg_method= x$bkg_method)
    }
    res$convergence <- ifelse(res$convergence==1, "convergence",
                        ifelse(res$convergence==2, "fit but no convergence",
                        ifelse(res$convergence==3, "error, no fit model",
                        NA)))
    return(res)
}

