##' @export
pars.survreg <- function(x,...) {
    c(coef(x),scale=x$scale)    
}


##' @export
score.survreg <- function(x,p,scale=TRUE,logscale=FALSE,indiv.logLik=FALSE,...) {    
    npar <- NROW(x$var)
    m <- model.frame(x)
    X <- model.matrix(terms(x), m)
    hasscale <- npar>length(x$coefficients)
    if (!missing(p)) {
        if (hasscale) sigma <- tail(p,1)
        p <- p[seq(length(p)-1)]
        x$linear.predictors <- as.vector(X%*%p)
        x$coefficients <- p
        x$scale <- sigma
    }
    derivatives <- residuals(x, type = "matrix")
    w <- model.weights(m)
    if (is.null(w)) w <- 1
    dldLP <- w*derivatives[,"dg"] ## Derivative wrt linear-predictor p=Xbeta
    S <- apply(X,2,function(x) x*dldLP)
    if (!is.null(x$naive.var)) {
        V <- x$naive.var
    } else {
        V <- x$var
    }
    if (hasscale && scale) {
        ds <- cbind("logsigma"=derivatives[,"ds"])
        if (!logscale) {
            ds <- ds/x$scale
            names(ds) <- "sigma"            
        }
        S <- cbind(S,ds)
    }
    if (hasscale && !scale) {
        V <- V[-npar,-npar,drop=FALSE]
    }
    attributes(S)$logLik <- 
                    if (indiv.logLik) derivatives[,"g"]
                    else sum(derivatives[,"g"])    
    attributes(S)$bread <- V
    return(S)
}

