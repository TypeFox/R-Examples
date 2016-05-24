#' Calculate cross-validated c-index
#' 
#' This function calculates cross-validated versions of Harrell's c-index.
#' 
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param type One of \code{"single"}, \code{"pair"} or \code{"fullpairs"}. For
#' \code{"single"} (default), the prognostic index Z_i is replaced by Z_i,(-i),
#' for \code{"pair"}, two assessments of concordance are made for each pair
#' (i,j), one using Z_i,(-i) and Z_j,(-i), the other using Z_i,(-j) and
#' Z_j,(-j), for \code{"fullpairs"}, each of the possible pairs is left out and
#' comparison is based on Z_i,(-i,-j) and Z_j,(-i,-j)
#' @param matrix if \code{TRUE}, the matrix of cross-validated prognostic
#' indices is also returned; default is \code{FALSE}
#' @return A list with elements \item{concordant}{The number of concordant
#' pairs} \item{total}{The total number of pairs that can be evaluated}
#' \item{cindex}{The cross-validated c-index} \item{matrix}{Matrix of
#' cross-validated prognostic indices (only if argument \code{matrix} is
#' \code{TRUE}} and with attribute \code{"type"} as given as input.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Harrell FE, Lee KL & Mark DB (1996), Multivariable prognostic
#' models: issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors, Statistics in Medicine 15, 361-387.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' # Real thing takes a long time, so on a smaller data set
#' ova2 <- ova[1:100,]
#' # Actual c-index
#' cindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2)
#' # Cross-validated c-indices
#' CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2)
#' CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2,
#'          type="pair")
#' \donttest{
#' CVcindex(Surv(tyears,d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova2,
#'          type="fullpairs")
#' }
#' 
#' @importFrom stats as.formula
#' 
#' @export CVcindex
CVcindex <- function(formula, data, type="single", matrix=FALSE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    # Center covariates
    # First get design matrix
    X <- model.matrix(formula, data = data)
    X <- X[,-1] # remove intercept
    X <- t(t(X) - apply(X,2,mean))
    cfull <- coxph(Surv(time,status) ~ X, data = data, method="breslow")
    n <- nrow(data)
    if (type=="single" | type=="pair") {
        m <- floor(log10(n))+1 # for formatting progress count
        pre <- rep("\b",2*m+1)
        xmat <- matrix(NA,n,n) # in columns (-i)
        # xmat[j,i] will contain PI_{j,(-i)}
        for (i in 1:n) {
            cat(pre,i,"/",n,sep=""); flush.console()
            # leave out i
            cmin1 <- coxph(Surv(time[-i], status[-i]) ~ X[-i,], method="breslow")
            # evaluate at all j except i
            xmat[-i,i] <- cmin1$linear.predictors
            # evaluate at i
            xmat[i,i] <- sum(X[i,] * cmin1$coef)
        }
        cat("\n")
    }
    if (type=="single") {
        formula <- as.formula("Surv(time,status) ~ x")
        ndata <- data.frame(time=time,status=status,x=diag(xmat))
        res <- cindex(formula=formula, data=ndata)
        if (matrix) res <- list(concordant=res$concordant,total=res$total,cindex=res$cindex,matrix=xmat)
    }
    if (type=="pair") {
        n <- length(time) # check if = length(status) and length(x)
        ord <- order(time,-status)
        time <- time[ord]
        status <- status[ord]
        xmat <- xmat[ord,ord]
        # pairs (i,j) for which the smallest observed time is an event time
        wh <- which(status==1)
        total <- concordant <- 0
        for (i in wh) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    if (time[j] > time[i]) { # ties not counted
                        total <- total + 2
                        if (xmat[j,i] < xmat[i,i]) concordant <- concordant + 1
                        if (xmat[j,j] < xmat[i,j]) concordant <- concordant + 1
                    }
                }
            }
        }
        if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
    }
    if (type=="fullpairs") {
        m <- floor(log10(n*(n-1)/2))+1 # for formatting progress count
        pre <- rep("\b",2*m+1)
        # xmat[i,j] will contain PI_{i,(-i,-j)}; xmat[j,i] will contain PI_{j,(-i,-j)}
        xmat <- matrix(NA,n,n)
        cnt <- 0
        for (i in 1:n) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    cnt <- cnt+1
                    cat(pre,cnt,"/",n*(n-1)/2,sep=""); flush.console()
                    # leave out i and j
                    cmin2 <- coxph(Surv(time[-c(i,j)], status[-c(i,j)]) ~ X[-c(i,j),], method="breslow")
                    # evaluate at i
                    xmat[i,j] <- sum(X[i,] * cmin2$coef)
                    # evaluate at j
                    xmat[j,i] <- sum(X[j,] * cmin2$coef)
                }
            }
        }
        ord <- order(time,-status)
        time <- time[ord]
        status <- status[ord]
        xmat <- xmat[ord,ord]
        # pairs (i,j) for which the smallest observed time is an event time
        wh <- which(status==1)
        total <- concordant <- 0
        for (i in wh) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    if (time[j] > time[i]) {# ties not counted
                        total <- total + 1
                        if (xmat[j,i] < xmat[i,j]) concordant <- concordant + 1
                    }
                }
            }
        }
        if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
    }
    cat("\n")
    attr(res, "type") <- type
    return(res)
}
