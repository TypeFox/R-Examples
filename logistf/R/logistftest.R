logistftest <-
function(object, test, values, firth = TRUE, beta0, weights, control)
{
   call <- match.call()
   formula<-object$formula
   data<-object$data
    if (missing(control)) control<-logistf.control()
    mf<-model.frame(object$formula,data=object$data)
##    mf <- match.call(expand.dots =FALSE)
 ##   m <- match(c("formula", "data","weights", "na.action", 
 ##       "offset"), names(mf), 0L)
## #   mf<-model.frame(formula, data=data, weights=weights)
##    mf <- mf[c(1, m)]
##    mf$drop.unused.levels <- TRUE
##    mf[[1L]] <- as.name("model.frame")
##    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    cov.name <- labels(x)[[2]]
#    weight <- as.vector(model.weights(mf)  )
    if(missing(weights) & !is.null(object$weights)) weight <- object$weights
    else weight<-NULL
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    
    cov.name <- labels(x)[[2]]
    k <- ncol(x)
    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }

###    fit.full<-logistf.fit(    ) # unrestricted, define init and col.fit from values, beta0 and test
###    fit.null<-logistf.fit(    ) # restricted, define init and col.fit from values, beta0 and test
    
    fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, control=control)
    
    pos<-coltotest
    if(missing(test))
        test <- coltotest
    if(is.vector(test))
        cov.name2 <- cov.name[test]
    else cov.name2 <- labels(model.matrix(test, data = data))[[2]]
    pos <- match(cov.name2, cov.name)   ## Position der Testfakt.
    OK <- !is.na(pos)
    pos <- pos[OK]
    cov.name2 <- cov.name2[OK]
    k2 <- length(cov.name2) ## Anzahl Faktoren
    if(!missing(beta0))
        offset1 <- beta0
    else offset1 <- rep(0, k)    ## Vektor der fixierten Werte
    if(!missing(values))
        offset1[pos] <- values
    beta <- offset1  ########################################

    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=(1:k)[-pos], control=control, init=beta)

    loglik<-c(fit.null$loglik,fit.full$loglik)
    
    offset1[ - pos] <- NA
    names(offset1) <- cov.name
    fit <- list(testcov = offset1, loglik = loglik, df = k2, prob = 1 - pchisq(2 *
        diff(loglik), k2), call = match.call(), beta = beta)
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    attr(fit, "class") <- "logistftest"
    fit
}

