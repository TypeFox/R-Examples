#AIC <- function(object, ...)
#    UseMethod("AIC")
"AIC.car" <-
function (object, ..., k=NULL) 
{
    cat("\nCall:\n", deparse(object$call), "\n", sep = "")
    phi <- object$phi
    covmat <- object$ecov
    ARP <- length(phi)
    tphi <- chol(solve(covmat)) %*% phi
    aic <- NULL
    cums <- cumsum(tphi^2)
    aic <- -cums + 2 * seq(1:length(phi))
    tphi <- round(as.vector(tphi), 2)
    aic <- round(aic, 2)
    res <- data.frame("order"=1:length(tphi), "t.statistic" = tphi, AIC = aic)
    cat("\nModel selection statistics", "\n\n")
    print(res, row.names=FALSE)
}
