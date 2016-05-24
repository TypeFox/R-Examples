
###{{{ distribution

##' @export
"distribution<-" <- function(x,...,value) UseMethod("distribution<-")

##' @export
"distribution" <- function(x,...,value) UseMethod("distribution")

##' @export
"distribution<-.lvm" <- function(x,variable,parname=NULL,init.par,mdist=FALSE,...,value) {
  if (inherits(variable,"formula")) variable <- all.vars(variable)
  dots <- list(...)
  ## Generator <- (is.function(value) && inherits(try(do.call(value,list()),silent=TRUE),"try-error"))
  ## if (Generator && length(variable)==1) value <- list(value)
  if (!is.null(parname) || length(dots)>0) { ## || Generator) {
      if (length(parname)>1 || (is.character(parname))) {
          if (missing(init.par)) {
              parameter(x,start=rep(1,length(parname))) <- parname
          } else {
              parameter(x,start=init.par) <- parname
          }
          ## if ("..."%ni%names(formals(value))) formals(value) <- c(formals(value),alist(...=))
          ## formals(value) <- modifyList(formals(value),dots)
          gen <- function(n,p,...) {
              args <- c(n,as.list(p[parname]),dots)
              names(args) <- names(formals(value))[seq(length(parname)+1)]
              do.call(value,args)
          }
      } else {
          gen <- value
          if ("..."%ni%names(formals(gen))) formals(gen) <- c(formals(gen),alist(...=))
          formals(gen) <- modifyList(formals(gen),dots)
          ## gen <- function(n,p,...) {
          ##     args <- c(n=n,dots)
          ##     names(args)[1] <- names(formals(value))[1]
          ##     do.call(value,args)
          ## }
      }
      ##      if (length(variable)>1)
          {
          gen <- list(gen)
      }
      distribution(x,variable,mdist=TRUE) <- gen
      return(x)
  }

  if (length(variable)==1 && !mdist) {
      addvar(x) <- as.formula(paste("~",variable))
      if (is.numeric(value)) value <- list(value)
      if (!is.null(attributes(value)$mean)) intercept(x,variable) <- attributes(value)$mean
      if (!is.null(attributes(value)$variance)) variance(x,variable) <- attributes(value)$variance
##      if (is.function(value) && "..."%ni%names(formals(value))) formals(value) <- c(formals(value),alist(...=))
      x$attributes$distribution[[variable]] <- value
      ## Remove from 'mdistribution'
      vars <- which(names(x$attributes$mdistribution$var)%in%variable)
      for (i in vars) {
          pos <- x$attributes$mdistribution$var[[i]]
          x$attributes$mdistribution$fun[pos] <- NULL
          x$attributes$mdistribution$var[which(x$attributes$mdistribution$var==pos)] <- NULL
          above <- which(x$attributes$mdistribution$var>pos)
          if (length(above)>0)
              x$attributes$mdistribution$var[above] <- lapply(x$attributes$mdistribution$var[above],function(x) x-1)
      }
      return(x)
  }

  if (is.list(value) && length(value)==1 && (is.function(value[[1]]) || is.null(value[[1]]))) {
      addvar(x) <- variable
      ## Multivariate distribution
      if (is.null(x$attributes$mdistribution)) x$attributes$mdistribution <- list(var=list(), fun=list())
      vars <- x$attributes$mdistribution$var

      if (any(ii <- which(names(vars)%in%variable))) {
          num <- unique(unlist(vars[ii]))
          vars[which(unlist(vars)%in%num)] <- NULL
          newfunlist <- list()
          numleft <- unique(unlist(vars))
          for (i in seq_along(numleft)) {
              newfunlist <- c(newfunlist, x$attributes$mdistribution$fun[[numleft[i]]])
              ii <- which(unlist(vars)==numleft[i])
              vars[ii] <- i
          }
          K <- length(numleft)
          x$attributes$mdistribution$var <- vars
          x$attributes$mdistribution$fun <- newfunlist
      } else {
          K <- length(x$attributes$mdistribution$fun)
      }
      if (length(distribution(x))>0)
          distribution(x,variable) <- rep(list(NULL),length(variable))

      x$attributes$mdistribution$var[variable] <- K+1
      x$attributes$mdistribution$fun <- c(x$attributes$mdistribution$fun,value)

      return(x)
  }

  if ((length(value)!=length(variable) & length(value)!=1))
      stop("Wrong number of values")
##  if (length(value)==1 && "..."%ni%names(formals(value))) formals(value) <- c(formals(value),alist(...=))
  for (i in seq_along(variable))
      if (length(value)==1) {
      distribution(x,variable[i],...) <- value
  } else {
##      if ("..."%ni%names(formals(value[[i]]))) formals(value[[i]]) <- c(formals(value[[i]]),alist(...=))
      distribution(x,variable[i],...) <- value[[i]]
  }
  return(x)

}

##' @export
"distribution.lvm" <- function(x,var,value,multivariate=FALSE,...) {
    if (!missing(value)) {
        distribution(x,var,...) <- value
        return(x)
    }
    if (multivariate) return(x$attributes$mdistribution)
  x$attributes$distribution[var]
}

###}}} distribution

###{{{ normal/gaussian

##' @export
normal.lvm <- function(link="identity",mean,sd,log=FALSE,...) {
  rnormal <- if(log) rlnorm else rnorm
  fam <- stats::gaussian(link); fam$link <- link
  f <- function(n,mu,var,...) rnormal(n,fam$linkinv(mu),sqrt(var))
  if (!missing(mean)) attr(f,"mean") <- mean
  if (!missing(sd)) attr(f,"variance") <- sd^2
  attr(f,"family") <- fam
  return(f)
}

##' @export
gaussian.lvm <- normal.lvm

##' @export
lognormal.lvm <- function(...) normal.lvm(...,log=TRUE)


###}}} normal/gaussian

###{{{ poisson

##' @export
poisson.lvm <- function(link="log",lambda,...) {
    fam <- stats::poisson(link); fam$link <- link
    f <- function(n,mu,...) {
        if (missing(n)) {
            return(fam)
        }
        rpois(n,fam$linkinv(mu))
    }
    if (!missing(lambda)) attr(f,"mean") <- fam$linkfun(lambda)
    attr(f,"family") <- fam
    attr(f,"var") <- FALSE
    return(f)
}

###}}} poisson

###{{{ pareto

## @examples
## m <- lvm()
## categorical(m,K=3) <- ~x
## distribution(m,~y) <- pareto.lvm(lambda=1)
## regression(m,additive=FALSE) <- y~x
## regression(m) <- y~z
## d <- sim(m,1e4,p=c("y~x:0"=1,"y~x:1"=1,"y~x:2"=exp(1)))
##
## X <- model.matrix(y~-1+factor(x)+z,data=d)
## mlogL <- function(theta) {
##     lambda <- exp(theta[1])
##     mu <- exp(X%*%theta[-1])
##     -sum(log(lambda*mu*(1+mu*d$y)^{-lambda-1}))
## }
## nlminb(rep(0,ncol(X)+1),mlogL)
##' @export
pareto.lvm <- function(lambda=1,...) {   ## shape: lambda, scale: mu
    ## Density f(y): lambda*mu*(1+mu*y)^{-lambda-1}
    ## Survival S(y): (1+mu*y)^{-lambda}
    ## Inverse CDF: u -> ((1-u)^{-1/lambda}-1)/mu
    f <- function(n,mu,var,...) {
        ((1-runif(n))^(-1/lambda)-1)/exp(mu)
    }
    attr(f,"family") <- list(family="pareto",
                             par=c(lambda=lambda))
    return(f)
}

###}}} poisson

###{{{ threshold

##' @export
threshold.lvm <- function(p,labels=NULL,...) {
    if (sum(p)>1 || any(p<0 | p>1)) stop("wrong probability vector") ;
    if (!is.null(labels))
    return(function(n,...) {
        return(cut(rnorm(n),breaks=c(-Inf,qnorm(cumsum(p)),Inf),labels=labels))
    })
    function(n,...)
        cut(rnorm(n),breaks=c(-Inf,qnorm(cumsum(p)),Inf))
}

###}}}

###{{{ binomial

##' @export
binomial.lvm <- function(link="logit",p,size=1) {
    fam <- stats::binomial(link); fam$link <- link
    f <- function(n,mu,var,...) {
        if (missing(n)) {
            return(fam)
        }
        rbinom(n,size,fam$linkinv(mu))
    }
    attr(f,"family") <- fam
    attr(f,"var") <- FALSE
    if (!missing(p)) attr(f,"mean") <- fam$linkfun(p)
    ## f <- switch(link,
    ##             logit =
    ##             function(n,mu,var,...) rbinom(n,1,tigol(mu)),
    ##             cloglog =
    ##             function(n,mu,var,...) rbinom(n,1,1-exp(-exp(1-mu))),
    ##             function(n,mu,var=1,...) rbinom(n,1,pnorm(mu,sd=sqrt(var)))
    ##             ### function(n,mu=0,var=1,...) (rnorm(n,mu,sqrt(var))>0)*1
    ##             )
    ##}
    return(f)
}

##' @export
logit.lvm <- binomial.lvm("logit")

##' @export
probit.lvm <- binomial.lvm("probit")

###}}} binomial

###{{{ Gamma


##' @export
Gamma.lvm <- function(link="inverse",shape,rate,unit=FALSE,var=FALSE,log=FALSE,...) {
  fam <- stats::Gamma(link); fam$link <- link
  rgam <- if (!log) rgamma else function(...) log(rgamma(...))
  if (!missing(shape) & !missing(rate))
    f <- function(n,mu,var,...) rgam(n,shape=shape,rate=rate)
  if (!missing(shape) & missing(rate)) {
    if (unit)
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=shape)
    else if (var)
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=sqrt(shape/var))
    else
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=shape/fam$linkinv(mu))
  }
  if (missing(shape) & !missing(rate)) {
    if (unit)
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=rate)
    else if (var)
      f <- function(n,mu,var,...) rgam(n,shape=rate^2*var,rate=rate)
    else
      f <- function(n,mu,var,...) rgam(n,shape=rate*fam$linkinv(mu),rate=rate)
  }
  if (missing(shape) & missing(rate)) {
    if (var)
      f <- function(n,mu,var,...) rgam(n,shape=var,rate=1)
    else
      f <- function(n,mu,var,...) rgam(n,shape=fam$linkinv(mu),rate=1)
  }
  attr(f,"family") <- fam
  attr(f,"var") <- FALSE
  return(f)
}

##' @export
loggamma.lvm <- function(...) Gamma.lvm(...,log=TRUE)


###}}} Gamma

###{{{ chisq
##' @export
chisq.lvm <- function(df=1,...) {
    function(n,mu,var,...) mu + rchisq(n,df=df)
}
###}}} chisq

###{{{ student (t-distribution)

##' @export
student.lvm <- function(df=2,mu,sigma,...) {
    f <- function(n,mu,var,...) mu + sqrt(var)*rt(n,df=df)
    if (!missing(mu)) attr(f,"mean") <- mu
    if (!missing(sigma)) attr(f,"variace") <- sigma^2
    return(f)
}

###}}} student (t-distribution)

###{{{ uniform

##' @export
uniform.lvm <- function(a,b) {
  if (!missing(a) & !missing(b))
    f <- function(n,mu,var,...) mu+runif(n,a,b)
  else
    f <- function(n,mu,var,...)
      (mu+(runif(n,-1,1)*sqrt(12)/2*sqrt(var)))
  return(f)
}

###}}} uniform

###{{{ weibull
## see also eventTime.R for coxWeibull

##' @export
weibull.lvm <- function(scale=100,shape=2) {
    ## accelerated failure time (AFT) regression
    ## parametrization.
    ##
    ## We parametrize the Weibull distribution (without covariates) as follows:
    ## hazard(t) = 1/shape * exp(-scale/shape) * t^(1/shape-1)
    ## The hazard is:
    ## - rising if shape > 1
    ## - declining if shape <1
    ## - constant if shape=1
    ##
    ## AFT regression
    ## hazard(t|Z) = 1/shape * exp(-scale/shape) * t^(1/shape-1) exp(-beta/shape*Z)
    ## scale^(-1/shape) = exp(a0+a1*X)
    ## PH regression
    ## scale = exp(b0+ b1*X)
    f <- function(n,mu,var,...) {
        ## message("weibull.scale=",scale)
        ## message("weibull.shape=",shape)
        ## message("coxWeibull.scale=",exp(scale/shape))
        (- log(runif(n)) * exp(scale/shape) * exp(mu/shape))^{shape}
        ## scale * (-log(1-runif(n)))^{1/shape}
        ## (- (log(runif(n)) / (1/scale)^(shape) * exp(-mu)))^(1/shape)
    }
    attr(f,"family") <- list(family="weibull",
                             regression="AFT",
                             par=c(shape=shape,scale=scale))
    return(f)
}

###}}} weibull

###{{{ sequence
##' @export
sequence.lvm <- function(a=0,b=1,integer=FALSE) {
    if (integer) {
        f <- function(n,...) seq(n)
        return(f)
    }
    if (is.null(a) || is.null(b)) {
        if (!is.null(a)) {
            f <- function(n,...) seq(a,length.out=n)
        } else {
            f <- function(n,...) seq(n)-(n-b)
        }
    } else {
        f <- function(n,...) seq(a,b,length.out=n)
    }
    return(f)
}
###}}} sequence

###{{{ ones
##' @export
ones.lvm <- function(p=1,interval=NULL) {
    f <- function(n,...) {
        if (!is.null(interval)) {
            val <- rep(0L,n)
            if (!is.list(interval)) interval <- list(interval)
            for (i in seq_along(interval)) {
                ii <- interval[[i]]
                lo <- round(ii[1]*n)
                hi <- round(ii[2]*n)
                val[seq(lo,hi)] <- 1L
            }
            return(val)
        }
        if (p==0) return(rep(0L,n))
        val <- rep(1L,n)
        if (p>0 && p<1) val[seq(n*(1-p))] <- 0L
        val
        }
  return(f)
}
###}}} ones

###}}}
