#' @method print gof.glm
#' @export
print.gof.glm <- function(x, ...){
    chiSq <- df <- pVal <- val <- test <- stat <- NULL
    m1 <- data.frame(x$chiSq[, list(chiSq, df, pVal)])
    rownames(m1) <- x$chiSq[, test]
    m1[, 1:3] <- format.pval(m1[, 1:3])
    stats::printCoefmat(m1,
                        has.Pvalue=TRUE,
                        cs.ind=1L,
                        dig.tst=getOption("digits"))
###
    m1 <- data.frame(x$gof[, list(val, df, pVal)])
    rownames(m1) <- x$gof[, paste(test, stat)]
    m1[is.na(m1)] <- NaN
    m1[, 1:3] <- format.pval(m1[, 1:3])
    suppressWarnings(
        stats::printCoefmat(m1,
                            has.Pvalue=TRUE,
                            cs.ind=1L,
                            dig.tst=getOption("digits")))
}
