#' @method print sig.glm
#' @export
print.sig.glm <- function(x, ...){
    pVal <- coef <- NULL
    m1 <- data.frame(
        "Wald"= x$Wald[, pVal],
        "LR"=x$LR[, pVal],
        "score"=x$score[, pVal])
    rownames(m1) <- x$Wald[, coef]
    m1[, 1:3] <- format.pval(m1[, 1:3])
    stats::printCoefmat(m1, has.Pvalue = T,
                        cs.ind=1L,
                        dig.tst=getOption("digits"))
}
