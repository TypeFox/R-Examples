##' Wrapper for automated backward elemination for logistic regression
##'
##' First run backward elemination via fastbw from the rms package, then
##' fit the logistic regression model
##' including the selected variables
##' @title Automated backward elemination for logistic regression
##' @param formula passed to lrm
##' @param data passed to lrm
##' @param ... passed to fastbw
##' @return object of class AutoSelectLRM
##' @seealso fastbw lrm
##' @examples
##' library(rms)
##' set.seed(7)
##' x <- abs(rnorm(20))
##' d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
##' fbw <- AutoSelectLRM(y~x+z,d)
##' predictStatusProb(fbw,newdata=d[1:3,])
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
AutoSelectLRM <- function(formula,data,...){
    fit <- rms::lrm(formula,data)
    bwfit <- rms::fastbw(fit,...)
    if (length(bwfit$names.kept)==0){
        newform <- reformulate("1",formula[[2]])
        newfit <- glm(newform,data,family="binomial")
    }
    else{
        newform <- reformulate(bwfit$names.kept, formula[[2]])
        newfit <- rms::lrm(newform,data,x=TRUE,y=TRUE)
    }
    out <- list(fit=newfit,In=bwfit$names.kept)
    out$call <- match.call()
    class(out) <- "AutoSelectLRM"
    out
}

#'@S3method predictStatusProb AutoSelectLRM
predictStatusProb.AutoSelectLRM <- function(object,newdata,...){
    p <- predictStatusProb(object[[1]],newdata=newdata,...)
    class(p) <- "predictStatusProb"
    p
}
