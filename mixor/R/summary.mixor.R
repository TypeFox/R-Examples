summary.mixor <-
function(object, digits=max(3, getOption("digits") - 2), signif.stars=TRUE, dig.tst = max(1, min(5, digits - 1)), ...) {
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
cat("Deviance =        ", object$Deviance, "\n")
cat("Log-likelihood = ", object$RLOGL,"\n")
cat("RIDGEMAX =        ", object$RIDGEMAX,"\n")
cat("AIC =            ", object$AIC, "\n")
cat("SBC =            ", object$SBC, "\n\n")
printCoefmat(object$Model, digits=digits, signif.stars=signif.stars, dig.tst=dig.tst, cs.ind=1:2, tst.ind=3, Pvalues=TRUE, has.Pvalue=TRUE)
}
