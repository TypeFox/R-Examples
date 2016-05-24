##' Wrapper function for glmnet 
##'
##' This function first calls \code{cv.glmnet} and then evaluates glmnet at the hyper parameter which optimizes the cross-validation criterion. 
##' @title Wrapper function for glmnet 
##' @param formula Formula where the right hand side specifies the
##' response and the left hand side the predictor matrix
##' @param data A data frame in which \code{formula} is evaluated
##' @param nfolds nfolds: number of cross-validation folds in
##' cv.glmnet (default in function is 10)
##' @param ... passed on to glmnet
##' @return
##' Object with class ElasticNet
##' @seealso predictStatusProb
##' @examples
##' # Generate some data with binary response Y
#'  # depending on X1 and X2 and X1*X2
#' set.seed(40)
#' N <- 40
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,1,.4)
#' X3 <- rnorm(N)
#' expit <- function(x) exp(x)/(1+exp(x))
#' lp <- expit(1 + X1 + X2 + X3)
#' Y <- factor(rbinom(N,1,lp))
#' dat <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
#'
#' efit <- ElasticNet(Y~X1+X2+X3,data=dat,family="binomial",alpha=0.1)
#' Brier(efit,verbose=FALSE)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
ElasticNet <- function(formula,
                       data,
                       nfolds=10,
                       ...){
    require(glmnet)
    call <- match.call()
    # get response and predictor variables
    mf <- model.frame(formula,data,na.action=na.omit)
    Terms <- terms.formula(formula)
    attr(Terms, "intercept") <- 1
    response <- model.response(mf)
    stopifnot(length(unique(response))==2)
    if (is.factor(response))
        response <- as.numeric(response==levels(response)[2])
    covariates <- model.matrix(Terms,data=mf)[,-1,drop=FALSE]
    # find lambda via cross-validation
    pathLambda <- glmnet::cv.glmnet(x=covariates,
                            y=response,
                            nfolds=nfolds,...)
    optlambda <- pathLambda$lambda.min
    # fit elastic net
    fit.enet <- glmnet::glmnet(x=covariates,
                       y=response,
                       lambda=optlambda,
                       ...)
    out <- list("call"=call, "enet"=fit.enet, "Lambda"=optlambda)
    class(out) <- "ElasticNet"
    invisible(out)
}

##' @S3method predictStatusProb ElasticNet
predictStatusProb.ElasticNet <- function(object, newdata,...){
    stopifnot(!missing(newdata))
    formula <- object$call$formula
    mf <- model.frame(formula,data=newdata,na.action=na.omit)
    terms <- attr(mf,"terms")
    newx <- model.matrix(terms,data=mf)[,-1]
    p <- predict(object$enet,
                 newx=newx,
                 type="response",
                 s=object$Lambda)
    as.numeric(p[,1,drop=TRUE])
    class(p) <- "predictStatusProb"
    p
}
