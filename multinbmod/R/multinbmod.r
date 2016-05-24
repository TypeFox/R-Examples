# Multivariate negative binomial model with robust estimation of regression coefficients.
# This function fits a multivariate negative binomial model by Maximum Likelihood and calculates robust standard errors
# of the regression coefficients.

multinbmod<-
function (formula, data, id, offset, start.coef = NULL,
start.phi = NULL,control=list())
# initial settings for optimisation
{   inv <- match.call()
    if (missing(data))
    data <- environment(formula)
       if (!is.list(control))
    stop("control must be a list")
    mf <- match.call(expand.dots = FALSE)
    mf$start.coef <- mf$start.phi<- mf$control<-NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, environment(formula))
    mt <- attr(mf, "terms")
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    p <- NCOL(X)
    Y <- model.response(mf, "numeric")

    id <- mf$"(id)"
    no.id <- (missing(id) || is.null(id))
    if (no.id) {
        warning("id variable is missing!")
        return(NULL)
    }

    offset <- model.offset(mf)
    if(is.null(offset)){
    offset <- rep(1,length(Y))}
    if (!is.null(offset)&& length(offset) != NROW(Y))
        stop(paste("Number of offsets is", length(offset), ", should equal",
            NROW(Y), "(number of observations)"))

if (is.null(start.coef)) {
        start.coef <- numeric(p)
        start.coef[1] <- log(mean(Y + 0.5))
        }

    else {
        if (length(start.coef) != p)
            stop("beta.start has wrong length")
    }
    if (is.null(start.phi)) {
        start.phi <- 0.5
    }
    else {
        if (length(start.phi) != 1)
            stop("start.sigma has wrong length")
        if (start.phi<=0)
            stop("start.phi must be positive")
           }

            fit<-multinb.fit(y=Y, x=X, offset, id=id, start.par=c(start.coef,start.phi),control)
            res<-list()
            res$converged <- fit$"converged?"
            res$coefficients<-fit$"estimated regression coefficients"
            res$model.coef.se<- fit$"se from model"
            res$robust.coef.se<- fit$"robust se"
            res$robust.t.values<- fit$"robust t-values"
            res$mle.phi <- fit$"estimated phi"
            res$phi.se<- fit$"se(phi)"
            res$minus2.loglik<-fit$"-2 x loglikelihood"
            res$iterations<- fit$iterations
            res$call<- inv
            class(res) <- "multinbmod"
            return(res)
            }


#Summary of a multinbmod object
summary.multinbmod<-function(object,...)
    {
    TAB<-cbind(Estimate=object$coefficients,
               ModelSE = object$model.coef.se,
               RobustSE = object$robust.coef.se,
               Robust.t= object$robust.t.values)
               res<- list( call=object$call,converged=object$converged,coefficients=TAB,MLE_of_phi=object$mle.phi,SE_of_phi=object$phi.se,
    minus2.loglik=object$minus2.loglik,iterations=object$iterations)
    class(res) <- "summary.multinbmod"
    res
    }


#Multivariate negative binomial model with robust estimation of regression coefficients
#This function is called by "multinbmod", but it can also be called directly
multinb.fit<-
function(y, x, offset=1, id, start.par,control=list())
{
    if (!is.list(control))
    stop("control must be a list")
        np<-NCOL(x)+1
        if (missing(start.par)) {
        start.par <- numeric(np)
        start.par[1] <- log(mean(y + 0.5))
        start.par[np]<- 0.5
        }
        else{
        if (length(start.par) != np)
            stop("start.par has wrong length")
        if (start.par[np]<=0)
            stop("starting value of phi must be positive")}
        no.id <- (missing(id) || is.null(id))
     if (no.id) {
       warning("id variable is missing!")
       return(NULL)
    }

      if (offset!=1 && length(offset) != NROW(y))
        stop(paste("Number of offsets is", length(offset), ", should equal",
            NROW(y), "(number of observations)"))
           tapplys.fun<-function(x,y){
           res<-tapply(x,y,sum)
           return(res[order(unique(y))])}

          ydot<-tapplys.fun(y, id)
        mnbloglik <- function(param, y, x, id,ydot)
        {
                p <- length(param)
                beta <- param[1:(p - 1)]
                a <- param[p]
                mu <- offset*exp(x %*% beta)
                mudot <- tapplys.fun(mu, id)
                minlog<- -1 * sum(lgamma(a^(-1) + ydot) - lgamma(a^(-1)) - (a^(-1) +
                        ydot) * log(1 + mudot * a) + ydot * log(a) + ydot *
                        log(mudot)) - sum(tapplys.fun(y * log(mu), id) - log(
                        mudot) * ydot)
                return(minlog)
        }

        grr <- function(param, y, x,id,ydot)
        {
                p <- length(param)
                beta <- param[1:(p - 1)]
                a <- param[p]
                mu <- offset*exp(x %*% beta)
                mudot <- tapplys.fun(mu, id)
                derivbeta <- rep(0, p)
                 for(i in 1:(p - 1)) {
                        derivbeta[i] <- -1*sum(tapplys.fun(y*x[,i],id) - tapplys.fun(x[,i]*mu,id)*(1+ydot*a)/(1+mudot*a))
                                                            }
                theta<- 1/a
                derivbeta[p] <- (1/a^2)*sum(digamma(theta+ydot)-log(1+mudot/theta)+
                (theta+ydot)*mudot/(theta*(theta+mudot))-ydot/theta-digamma(theta))
                return(derivbeta)
        }

      hes<-function(param,y,x,id,ydot)
        {
                p<- length(param)
                beta <- param[1:(p - 1)]
                a <- param[p]
                mu <- offset*exp(x %*% beta)
                mudot <- tapplys.fun(mu, id)
                hesmat<-matrix(numeric(p^2),p,p)
                  for (j in (1:(p-1))){
                  for (i in (1:(p-1))){
                   hesmat[i,j]<- -1*sum(-1*tapplys.fun(mu*x[,i]*x[,j],id)*(1+a*ydot)/(1+a*mudot) + tapplys.fun(mu*x[,i],id)*tapplys.fun(mu*x[,j],id)*a*
                   (1+a*ydot)/(1+a*mudot)^2)
                   }}
                   for (i in (1:(p-1))){
                hesmat[i,p]<-hesmat[p,i]<- sum(tapplys.fun(mu*x[,i],id)*(ydot-mudot)/(1+mudot*a)^2)
                }
                hesmat[p,p]<- -1*sum(-ydot/a^2 + trigamma(a^(-1)+ydot)/a^4+2*digamma(a^(-1)+ydot)/a^3-trigamma(a^(-1))/a^4-
       2*digamma(a^(-1))/a^3 + 2*mudot/(a^2*(1+mudot*a))-2*log(1+mudot*a)/a^3+(a^(-1)+ydot)*mudot^2/(1+mudot*a)^2)

       return(hesmat)

        }

        res<- nlminb(start.par, mnbloglik,grr,lower = c(rep( - Inf, (np - 1)), 0.01),
              y = y, x = x, id = id,ydot=ydot,control=control)

        betamle <- res$par[1:(np - 1)]
        names(betamle)<-names(as.data.frame(x))
        amle <- as.numeric(res$par[np])
        muest <- offset*exp(x %*% betamle)
        mudotest <- tapplys.fun(muest, id)

        svec <- (y - muest)
        Nid <- length(unique(id))
        var0<- var1 <- matrix(numeric((np - 1)^2), (np - 1))
        for(i in 1:Nid) {
                eks <- x[(pick <- id == unique(id)[i]),  , drop = F]
                eks <- as.matrix(eks)
                mus <- muest[pick]
                ders  <- eks * mus
                if(length(mus) > 1) {
                        visinv <- diag(1/mus) - (amle/(1 + amle * sum(mus))) * rep(1, length(mus)) %*% t(rep(1, length(
                                mus)))
                        var0 <- var0 + t(ders) %*% visinv %*% ders
                        vsan<- svec[pick] %*% t(svec[pick])
                        var1 <- var1 + t(ders) %*% visinv %*% vsan %*% visinv %*% ders
                }
                else {
                        visinv <- 1/(mus * (1 + amle * mus))
                        var0 <- var0 + visinv * (t(ders) %*% ders)
                        var1 <- var1 + (t(ders) %*% ders) * (visinv^2 * (svec[
                                pick])^2)
                }
        }

        modcov<-solve(var0)
        robcovbetamle <- solve(var0) %*% var1 %*% solve(var0)
        robsebetamle <- sqrt(diag(robcovbetamle))
        names(robsebetamle)<-names(betamle)
        tvalues <- betamle/robsebetamle
        modsebemle<-  sqrt(diag(modcov))
        names(modsebemle)<-names(betamle)
        amle<- res$par[np]
        varesta<- 1/(-1*sum(-ydot/amle^2 + trigamma(amle^(-1)+ydot)/amle^4+2*digamma(amle^(-1)+ydot)/amle^3-trigamma(amle^(-1))/amle^4-
       2*digamma(amle^(-1))/amle^3 + 2*mudotest/(amle^2*(1+mudotest*amle))-2*log(1+mudotest*amle)/amle^3+(amle^(-1)+ydot)*mudotest^2/(1+mudotest*amle)^2) )



 myres <- list(betamle,modsebemle,robsebetamle, tvalues,modcov,robcovbetamle,
                amle, sqrt(varesta), 2*(res$
                objective+sum(lgamma(y+1))), ifelse(res$convergence==0,TRUE,FALSE),res$iterations)

        names(myres) <- c("estimated regression coefficients","se from model","robust se", "robust t-values","covariance of beta estimates from model",
                "robust covariance of beta estimates","estimated phi","se(phi)","-2 x loglikelihood","converged?",
                 "iterations")
        myres
}




