summary.zlm <-
function (object, printout = TRUE, ...) 
{
    betas = object$coefficients
    betas2 = object$coef2moments
    sds = sqrt(betas2 - betas^2)
    ests = cbind(betas, sds)
    gi = object$gprior.info
    gi.choice = gi$gtype
    if (gi$gtype == "hyper") {
        gi.choice = paste(gi$gtype, " (a=", 2 + signif(gi$hyper.parameter - 
            2, digits = 4), ")", sep = "")
    }
    gi.sd = -1
    gi.sdtext = ""
    if (length(gi$shrinkage.moments) > 1) {
        gi.sd = sqrt(gi$shrinkage.moments[[2]] - gi$shrinkage.moments[[1]]^2)
        gi.sdtext = paste(" St.Dev.:", round(gi.sd, 3))
    }
    rownames(ests) = c("(Intercept)", attr(object$terms, "term.labels"))
    colnames(ests) = c("Exp.Val.", "St.Dev.")
    cat("Coefficients\n")
    print(ests)
    cat("\n Log Marginal Likelihood:\n")
    cat(object$marg.lik)
    cat(paste("\n g-Prior:", gi.choice, "\n"))
    cat(paste("Shrinkage Factor", ifelse(gi$is.constant, ": ", 
        " Exp.Val: "), round(gi$shrinkage.moments[[1]], 3), gi.sdtext, 
        "\n", sep = ""))
    res = list()
    res$residuals <- object$residuals
    res$coefficients <- object$coefficients
    res$coef.sd <- sds
    res$gprior <- gi.choice
    res$E.shrinkage <- gi$shrinkage.moments[[1]]
    if (gi.sd > -1) {
        res$SD.shrinkage <- gi.sd
    }
    res$log.lik <- object$marg.lik
    return(invisible(res))
}
