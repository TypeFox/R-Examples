#' Calculate cross-validated log-partial likelihood (with shrinkage)
#' 
#' This function calculates the cross-validated log partial likelihood, with
#' shrinkage if requested.
#' 
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param progress if \code{TRUE} (default), progress of the cross-validation
#' will be printed
#' @param overall if \code{TRUE}, \code{CVPL} uses regression coefficient
#' estimates based on the full data, for each observation i, rather than the
#' estimates based on data minus i
#' @param shrink Shrinkage factor; default is 1 (no shrinkage)
#' @return Numeric; the cross-validated log partial likelihood
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Verweij PJM & van Houwelingen HC (1994), Penalized likelihood in
#' Cox regression, Statistics in Medicine 13, 2427-2436.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' CVPL(Surv(tyears, d) ~ 1, data = ova)
#' CVPL(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam,
#'   data = ova)
#' CVPL(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam,
#'   data = ova, overall=TRUE)
#' 
#' @importFrom stats model.matrix
#' @importFrom utils flush.console
#' 
#' @export CVPL
CVPL <- function(formula, data, progress=TRUE, overall=FALSE, shrink=1)
{
    # Extract data (time and status)
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[, p - 1]
    status <- y[, p]
    n <- length(time)
    if (nrow(data) != n)
        stop("missing data in time or status not allowed")
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    # Model matrix
    X <- model.matrix(formula, data = data)
    if (ncol(X)>1) {
        X <- X[, -1, drop = FALSE]
        X <- X[ord, , drop = FALSE]
        X <- t(t(X) - apply(X,2,mean)) # Centering
    }
    if (overall) {
        cfull <- coxph(Surv(time,status) ~ X, data = data, method="breslow")
        pl1 <- cfull$loglik[2]
    }
    if (progress) {
        m <- floor(log10(n)) + 1
        pre <- rep("\b", 2 * m + 1)
    }
    rat <- rep(NA,n)
    res <- 0
    for (i in 1:n) {
        if (progress) {cat(pre,i,"/",n,sep=""); flush.console()}
        # leave out i
        if (!overall) { # estimate coefficients without i
            cmin1 <- coxph(Surv(time[-i], status[-i]) ~ X[-i,], method="breslow")
            rat[-i] <- exp(shrink * cmin1$linear.predictors)
            rat[i] <- exp(sum(X[i,] * shrink * cmin1$coef))
            rcsrat <- rev(cumsum(rev(rat)))
            rat1 <- rat[status==1]
            rcsrat1 <- rcsrat[status==1]
            pl1 <- sum(log(rat1) - log(rcsrat1))
            if (shrink==1) pl2 <- cmin1$loglik[2]
        }
        else # use overall estimate to define rat
            rat <- exp(cfull$linear.predictors)
        if (overall | shrink != 1) {
            # Can't use loglik of cmin1 now
            ratmin1 <- rat[-i]
            rcsratmin1 <- rev(cumsum(rev(ratmin1)))
            rat1 <- ratmin1[status[-i]==1]
            rcsrat1 <- rcsratmin1[status[-i]==1]
            pl2 <- sum(log(rat1) - log(rcsrat1))
        }
        res <- res + (pl1-pl2)
    }
    if (progress) cat(paste(rep("\b", 4*m+3), collapse=""))
    return(res)
}
