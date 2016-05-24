##' S3 wrapper function for party's ctree method
##'
##' The ModelGood crossvalidation functionality works only
##' for S3 classes. 
##' @title S3 wrapper function for party's ctree method
##' @param ... passed to \code{ctree}
##' @return object of class Ctree which contains a ctree object
##' @seealso ctree
##' @examples
##' library(party)
##' set.seed(7)
##' x <- abs(rnorm(20))
##' d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
##' ct <- Ctree(y~x+z,d)
##' plot(ct$ctree)
##' predictStatusProb(ct,newdata=d[1:3,])
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
Ctree <- function(...){
    ## require(party)
    out <- list(ctree=party::ctree(...))
    class(out) <- "Ctree"
    out$call <- match.call()
    out  
}

#'@S3method predictStatusProb Ctree
predictStatusProb.Ctree <- function (object, newdata, ...) {
    ## require(party)
    N <- NROW(newdata)
    p <- party::treeresponse(object$ctree, newdata=newdata)
    p <- sapply(p,function(x)x[1])
    if (NROW(p) != NROW(newdata))
        stop("Prediction failed")
    p
    class(p) <- "predictStatusProb"
    p
}
