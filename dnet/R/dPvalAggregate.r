#' Function to aggregate p values
#'
#' \code{dPvalAggregate} is supposed to aggregate a input matrix p-values into a vector of aggregated p-values. The aggregate operation is applied to each row of input matrix, each resulting in an aggregated p-value. The method implemented can be based on the order statistics of p-values or according to Fisher's method or Z-transform method.
#'
#' @param pmatrix a data frame or matrix of p-values
#' @param method the method used. It can be either "orderStatistic" for the method based on the order statistics of p-values, or "fishers" for Fisher's method, or "Ztransform" for Z-transform test (and the weighted Z-test)
#' @param order an integeter specifying the order used for the aggregation according to the order statistics of p-values
#' @param weight a vector specifying the weights used for the aggregation according to Z-transform method
#' @return 
#' \itemize{
#' \item{\code{ap}: a vector with the length nrow(pmatrix), containing aggregated p-values}
#' }
#' @note For each row of input matrix with the \eqn{c} columns, there are \eqn{c} p-values that are uniformly independently distributed over [0,1] under the null hypothesis (uniform distribution). According to the order statisitcs, they follow the Beta distribution with the paramters \eqn{a=order} and \eqn{b=c-order+1}. According to the Fisher's method, after transformation by \eqn{-2*\sum^clog(pvalue)}, they follow Chi-Squared distribution. According to the Z-transform method, first converts the one-tailed P-values into standard normal deviates Z, then combines Z via \eqn{\frac{\sum^c(w*Z)}{\sum^c(w^2)}}, where \eqn{w} is the weight (usually square root of the sample size if the weighted Z-test; 1 if Z-transform test), and finally the combined Z follows the standard normal distribution to test the cumulative/aggregated evidence on the common null hypothesis.
#' @export
#' @seealso \code{\link{dPvalAggregate}}
#' @include dPvalAggregate.r
#' @examples
#' # 1) generate an iid uniformly-distributed random matrix of 1000x3 
#' pmatrix <- cbind(runif(1000), runif(1000), runif(1000))
#'
#' # 2) aggregate according to the order statistics
#' ap <- dPvalAggregate(pmatrix, method="orderStatistic")
#'
#' # 3) aggregate according to the Fisher's method
#' ap <- dPvalAggregate(pmatrix, method="fishers") 
#'
#' # 4) aggregate according to the Z-transform method
#' ap <- dPvalAggregate(pmatrix, method="Ztransform")

dPvalAggregate <- function (pmatrix, method=c("orderStatistic", "fishers", "Ztransform"), order=ncol(pmatrix), weight=rep(1, ncol(pmatrix))) 
{
    
    method <- match.arg(method)
    
    if (is.vector(pmatrix)){
        pmatrix <- matrix(pmatrix, nrow=length(pmatrix), ncol=1)
    }else if(is.matrix(pmatrix) | is.data.frame(pmatrix)){
        pmatrix <- as.matrix(pmatrix)
    }else if(is.null(pmatrix)){
        stop("The input pmatrix must be not NULL.\n")
    }
    if(is.null(rownames(pmatrix))){
        rownames(pmatrix) <- seq(1,nrow(pmatrix))
    }
    nr <- nrow(pmatrix)
    nc <- ncol(pmatrix)
    
    if(method == "orderStatistic"){
        if(is.null(order)){
            order <- nc
        }
        order <- as.integer(order)
        if(order < 1){
            order <- 1
        }else if(order > nc){
            order <- nc
        }
        
        omatrix <- pmatrix
        for (j in 1:nr) {
            omatrix[j,] <- sort(as.numeric(pmatrix[j,]), na.last=T)
        }
        
        x <- as.numeric(omatrix[,order])
        ap <- stats::pbeta(x, shape1=order, shape2=nc-order+1) # distribution function for Beta distribution
        names(ap) <- rownames(pmatrix)
        
    }else if(method == "fishers"){
        
        fmatrix <- pmatrix
        # replace those zeros with the minumum of non-zeros
        fmatrix[fmatrix==0] <- min(fmatrix[fmatrix!=0])
        fmatrix <- log(fmatrix)
        x <- -2*apply(fmatrix, 1, sum)
        ap <- stats::pchisq(x, df=2*nc, lower.tail=F) # distribution function for Chi-Squared distribution
        names(ap) <- rownames(pmatrix)
        
    }else if(method == "Ztransform"){
        
        if(is.null(weight)){
            weight <- rep(1, ncol(pmatrix))
        }
        
        zmatrix <- (stats::qnorm(pmatrix, lower.tail=F) %*% weight)/sqrt(sum(weight^2))
        ap <- stats::pnorm(zmatrix, lower.tail=F) # distribution function for normal distribution
        names(ap) <- rownames(pmatrix)
        
    }

    return(ap)
}