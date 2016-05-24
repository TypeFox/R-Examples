#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: models.R 4909 2013-08-15 09:46:51Z gruen $
#

FLXMRglm <- function(formula=.~.,
                     family=c("gaussian", "binomial", "poisson", "Gamma"),
                     offset=NULL)
{
    family <- match.arg(family)
    glmrefit <- function(x, y, w) {
      fit <- c(glm.fit(x, y, weights=w, offset=offset,
                       family=get(family, mode="function")()),
               list(call = sys.call(), offset = offset,
                    control = eval(formals(glm.fit)$control),            
                    method = "weighted.glm.fit"))
      fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
      fit$df.residual <- sum(w) - fit$rank
      fit$x <- x
      fit
    }
                
    z <- new("FLXMRglm", weighted=TRUE, formula=formula,
             name=paste("FLXMRglm", family, sep=":"), offset = offset,
             family=family, refit=glmrefit)
    z@preproc.y <- function(x){
      if (ncol(x) > 1)
        stop(paste("for the", family, "family y must be univariate"))
      x
    }

    if(family=="gaussian"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <-  p + offset
          p
        }

        logLik <- function(x, y, ...)
          dnorm(y, mean=predict(x, ...), sd=sigma, log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=coef, sigma=sigma),
            logLik=logLik, predict=predict,
            df=df)
      })

      z@fit <- function(x, y, w, component){
        fit <- lm.wfit(x, y, w=w, offset=offset)
        with(list(coef = coef(fit), df = ncol(x)+1,
                  sigma =  sqrt(sum(fit$weights * fit$residuals^2 /
                    mean(fit$weights))/ (nrow(x)-fit$rank))),
             eval(z@defineComponent))
      }
    }
    else if(family=="binomial"){
      z@preproc.y <- function(x){
        if (ncol(x) != 2)
          stop("for the binomial family, y must be a 2 column matrix\n",
               "where col 1 is no. successes and col 2 is no. failures")
        if (any(x < 0))
          stop("negative values are not allowed for the binomial family")
        x
      }     
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dbinom(y[,1], size=rowSums(y), prob=predict(x, ...), log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=coef),
            logLik=logLik, predict=predict,
            df=df)
      })

      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=binomial(), offset=offset, start=component$coef)
        with(list(coef = coef(fit), df = ncol(x)),
             eval(z@defineComponent))
      }
    }
    else if(family=="poisson"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dpois(y, lambda=predict(x, ...), log=TRUE)
        
        new("FLXcomponent",
            parameters=list(coef=coef),
            logLik=logLik, predict=predict,
            df=df)
      })
          
      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=poisson(), offset=offset, start=component$coef)
        with(list(coef = coef(fit), df = ncol(x)),
             eval(z@defineComponent))
      }
    }
    else if(family=="Gamma"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dgamma(y, shape = shape, scale=predict(x, ...)/shape, log=TRUE)
        
        new("FLXcomponent", 
            parameters = list(coef = coef, shape = shape),
            predict = predict, logLik = logLik,
            df = df)
      })

      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=Gamma(), offset=offset, start=component$coef)
        with(list(coef = coef(fit), df = ncol(x)+1,
                  shape = sum(fit$prior.weights)/fit$deviance),
             eval(z@defineComponent))
      }
    }
    else stop(paste("Unknown family", family))
    z
}

###**********************************************************

FLXMCmvnorm <- function(formula=.~., diagonal=TRUE)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="model-based Gaussian clustering")

    z@defineComponent <- expression({
      logLik <- function(x, y)
        mvtnorm::dmvnorm(y, mean=center, sigma=cov, log=TRUE)
    
      predict <-  function(x, ...)
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      new("FLXcomponent", parameters=list(center = center, cov = cov),
          df=df, logLik=logLik, predict=predict)
    })
    
    z@fit <- function(x, y, w, ...){
      para <- cov.wt(y, wt=w)[c("center","cov")]
      para$df <- (3*ncol(y) + ncol(y)^2)/2
      if(diagonal){
        para$cov <- diag(diag(para$cov))
        para$df <- 2*ncol(y)
      }
      with(para,
           eval(z@defineComponent))
    }
    z
}

FLXMCnorm1 <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="model-based univariate Gaussian clustering")

    z@defineComponent <- expression({
      logLik <- function(x, y)
        dnorm(y, mean=center, sd=sqrt(cov), log=TRUE)
    
      predict <-  function(x, ...)
        matrix(center, nrow=nrow(x), ncol=1,
               byrow=TRUE)
      new("FLXcomponent",
          parameters=list(mean = as.vector(center), sd = as.vector(sqrt(cov))),
          df=df, logLik=logLik, predict=predict)
    })
    
    z@fit <- function(x, y, w, ...){
      para <- cov.wt(as.matrix(y), wt=w)[c("center","cov")]
      para$df <- 2
      with(para, eval(z@defineComponent))
    }
    z
}


###**********************************************************

FLXMCmvbinary <- function(formula=.~., truncated = FALSE) {
  if (truncated) return(MCmvbinary_truncated(formula))
  else return(MCmvbinary(formula))
}

MCmvbinary <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvbinary", name="model-based binary clustering")

    ## make sure that y is binary
    z@preproc.y <- function(x){
        storage.mode(x) <- "logical"
        storage.mode(x) <- "integer"
        x
    }
    z@defineComponent <- expression({
      predict <- function(x, ...){
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      }
        
      logLik <- function(x, y){
        p <- matrix(center, nrow=nrow(x), ncol=length(center),
                    byrow=TRUE)
        rowSums(log(y*p+(1-y)*(1-p)))
      }
            
      new("FLXcomponent", parameters=list(center=center), df=df,
          logLik=logLik, predict=predict)
    })

    z@fit <- function(x, y, w, ...)
      with(list(center = colSums(w*y)/sum(w), df = ncol(y)),
           eval(z@defineComponent))
    
    z
}




###**********************************************************

binary_truncated <- function(y, w, maxit = 200, epsilon = .Machine$double.eps) {
  r_k <- colSums(y*w)/sum(w)
  r_0 <- 0
  llh.old <- -Inf
  for (i in seq_len(maxit)) {
    p <- r_k/(1+r_0)
    llh <- sum((r_k*log(p))[r_k > 0])+ sum(((1 - r_k + r_0) * log(1-p))[(1-r_k+r_0) > 0])
    if (abs(llh - llh.old)/(abs(llh) + 0.1) < epsilon) break    
    llh.old <- llh
    prod_p <- prod(1-p)
    r_0 <- prod_p/(1-prod_p)
  }
  p
}

MCmvbinary_truncated <- function(formula=.~.)
{
    z <- MCmvbinary(formula=formula)
    z@defineComponent <- expression({
      predict <- function(x, ...) {
        matrix(center, nrow = nrow(x), ncol = length(center), 
               byrow = TRUE)
      }
      logLik <- function(x, y) {
        p <- matrix(center, nrow = nrow(x), ncol = length(center), 
                    byrow = TRUE)
        rowSums(log(y * p + (1 - y) * (1 - p))) - log(1 - prod(1-center))
      }
      new("FLXcomponent", parameters = list(center = center), df = df, 
          logLik = logLik, predict = predict)
    })
    z@fit <- function(x, y, w, ...){
      with(list(center = binary_truncated(y, w), df = ncol(y)),
           eval(z@defineComponent))
    }   
    z
}


###**********************************************************

FLXMCmvcombi <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvcombi",
             name="model-based binary-Gaussian clustering")

    ## figure out who is binary
    BINARY <- NULL
    z@preproc.y <- function(x){
      BINARY <<- apply(x, 2, function(z) all(unique(z) %in% c(0,1)))
      x
    }
    
    z@defineComponent <- expression({
      predict <- function(x, ...){
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      }
      
      logLik <- function(x, y){
        z <- 0
        if(any(BINARY)){
          p <- matrix(center[BINARY], nrow=nrow(x),
                      ncol=sum(BINARY), byrow=TRUE)
          z <- z+rowSums(log(y[,BINARY,drop=FALSE]*p +
                             (1-y[,BINARY,drop=FALSE])*(1-p)))
        }
        if(!all(BINARY)){
          if(sum(!BINARY)==1)
            z <- z + dnorm(y[,!BINARY],
                           mean=center[!BINARY], sd=sqrt(var[!BINARY]),
                           log=TRUE)
          else
            z <- z + mvtnorm::dmvnorm(y[,!BINARY,drop=FALSE],
                                      mean=center[!BINARY], sigma=diag(var[!BINARY]),
                                      log=TRUE)
        }
        z
      }
            
      new("FLXcomponent", parameters=list(center=center, var=var), df=df,
          logLik=logLik, predict=predict)
    })

    z@fit <- function(x, y, w, ...){
      para <- cov.wt(y, wt=w)[c("center","cov")]
      para$var <- diag(para$cov)
      para$var <- pmax(para$var, sqrt(.Machine$double.eps))
      para$var[BINARY] <- NA
      para$cov <- NULL
      para$df <- ncol(y) + sum(!BINARY)


      with(para, eval(z@defineComponent))
    }
    z
}

