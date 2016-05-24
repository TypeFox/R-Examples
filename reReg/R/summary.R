print.reReg <- function(x, ...) {
    if (!is.reReg(x))
        stop("Must be a reReg object")
    cat("Call:\n")
    print(x$Call)
    if (x$method == "jsc") 
        ## cat("\nMethod:", x$method, "\n")
        cat("\nMethod: Joint Scale-Change Model \n")
    if (x$method == "GL")
        cat("\nMethod: Ghosh-Lin method \n")
    if (x$method == "HW")
        cat("\nMethod: Huang-Wang method \n")
    cat("\nCoefficients:\n")
    coefMat <- rbind(x$alpha, x$beta)
    rownames(coefMat) <- c("alpha", "beta")
    print(coefMat)
}

summary.reReg <- function(object, ...) {
    if (!is.reReg(object))
        stop("Must be a reReg object")
    tabA <- cbind(Estimate = round(object$alpha, 3),
                  StdErr = round(object$aSE, 3),
                  z.value = round(object$alpha / object$aSE, 3),
                  p.value = round(2 * pnorm(-abs(object$alpha / object$aSE)), 3))
    rownames(tabA) <- names(object$alpha)
    tabB <- cbind(Estimate = round(object$beta, 3),
                  StdErr = round(object$bSE, 3),
                  z.value = round(object$beta / object$bSE, 3),
                  p.value = round(2 * pnorm(-abs(object$beta / object$bSE)), 3))
    rownames(tabB) <- names(object$beta)
    out <- list(Call = object$Call, method = object$method, tabA = tabA, tabB = tabB)
    class(out) <- "summary.reReg"
    out
}

print.summary.reReg <- function(x, ...) {
    cat("Call:\n")
    print(x$Call)
    if (x$method == "jsc") 
        cat("\nMethod: Joint Scale-Change Model \n")
    if (x$method == "GL")
        cat("\nMethod: Ghosh-Lin method \n")
    if (x$method == "HW")
        cat("\nMethod: Huang-Wang method \n")
    cat("\nCoefficients (rate):\n")
    printCoefmat(as.data.frame(x$tabA), P.values = TRUE, has.Pvalue = TRUE)
    cat("\nCoefficients (hazard):\n")
    printCoefmat(as.data.frame(x$tabB), P.values = TRUE, has.Pvalue = TRUE)
}
