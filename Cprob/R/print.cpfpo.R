print.cpfpo <- function(x,...) {
    if (!inherits(x,"cpfpo"))
        stop ("'x' must be of class 'cpfpo'")
    ncov <- ncol(x$alpha)
    namen <- colnames(x$alpha)
    cat("Test for non-significant effect \n")
    sig <- sig.test.int.ff(x,idx=1:ncov,ncut=0)
    rownames(sig) <- namen
    colnames(sig) <- c("stat","pvalue")
    print(sig)
    cat("\n Test for constant fit \n")
    const <- cst.fit.ff(x,1:ncov)$eta
    test <- gof.test.int.ff(x,idx=1:ncov,ncut=0)
    res <- data.frame(coef=const$eta,
                      exp.coef=exp(const$eta),
                      SE.coef=const$eta.se,
                      stat=test[,1],
                      pvalue=test[,2]
                      )
    rownames(res) <- namen
    print(res)
    invisible(list(sig = sig, tdep = res))
}
