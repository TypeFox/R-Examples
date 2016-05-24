##' The call is added to an ctree object
##' 
##' @title S3-Wrapper for ctree.
##' @param ... passed to ctree
##' @return list with two elements: ctree and call
##' @seealso pecCforest
##' @examples
##' library(prodlim)
##' library(party)
##' library(survival)
##' set.seed(50)
##' d <- SimSurv(50)
##' nd <- data.frame(X1=c(0,1,0),X2=c(-1,0,1))
##' f <- pecCtree(Surv(time,status)~X1+X2,data=d)
##' predictSurvProb(f,newdata=nd,times=c(3,8))
##' 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export 
pecCtree <- function(...){
 out <- list(ctree=party::ctree(...))
 class(out) <- "pecCtree"
 out$call <- match.call()
 out  
}

##' @export 
predictSurvProb.pecCtree <- function (object, newdata, times, ...) {
    requireNamespace("party")
    N <- NROW(newdata)
    NT <- length(times)
    survObj <- party::treeresponse(object$ctree, newdata=newdata)
    p <- do.call("rbind", lapply(survObj,function(x){
        predictSurvProb(x, newdata=newdata[1,,drop=FALSE], times=times)
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}
