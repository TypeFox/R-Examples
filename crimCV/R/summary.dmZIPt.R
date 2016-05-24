summary.dmZIPt <-
function (object, ...) 
{
    cat("Summary \n")
    cat("======= \n \n")
    cat("AIC: ", object$AIC, "\n")
    cat("BIC: ", object$BIC, "\n")
    if (!is.null(object$cv)) {
        cat("CVE: ", object$cv, "\n")
    }
    gwt <- object$gwt
    gwt <- round(gwt, 3)
    gwt <- gwt/apply(gwt, 1, sum)
    colnames(gwt) <- paste("Grp", 1:object$ng, sep = "")
    rownames(gwt) <- paste("Indiv", 1:object$ni, sep = "")
    cat("\n")
    cat("Group Membership Prob. \n")
    cat("====================== \n \n")
    print(gwt)
}
