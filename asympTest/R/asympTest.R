asymp.test <- function(x, ...) UseMethod("asymp.test")

asymp.test.default <- function(x, y = NULL,parameter=c("mean","var","dMean","dVar","rMean","rVar"),
                         alternative = c("two.sided", "less", "greater"),
                               reference=0, conf.level = 0.95, rho=1,...) {

  parameter <- match.arg(parameter)
  alternative <- match.arg(alternative)

  if(parameter=="mean" && !is.null(y)) parameter <- "dMean"

  if (parameter %in% c("dMean","dVar","rMean","rVar") & is.null(y)) stop(" 'y' is missing for two-sample test")


  if ((parameter %in% c("mean","var","dMean","dVar","rMean","rVar"))==FALSE) stop(" 'parameter' must be one of 'mean','var','dMean','dVar','rMean','rVar'")


  if (is.null(y)) data.name <- deparse(substitute(x))

  if (!is.null(y)){
    if (parameter%in%c("mean","var")) data.name <- deparse(substitute(x))
    else data.name <- paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="")
  }	
  
  
  if ((alternative %in% c("two.sided", "less", "greater"))==FALSE) stop(" 'alternative' must be one of 'two.sided', 'less', or 'greater'")
  
  if (!is.numeric(reference)) stop(" 'reference' must be numeric")
  
  
  if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1)) stop("'conf.level' must be a single number between 0 and 1")
  
  if (parameter %in%c("mean","var")){
    if (length(x)<2) stop(" 'x' must be at least of length >1")
  }
  else {
    if ((length(x)<2) | (length(y)<2)) stop(" 'x' and 'y' must be at least of length >1")
  }
  
  
  switch(parameter,
         mean={ 
           thetaEst <- mean(x)
           seThetaEst <- seMean(x)
           txtPar <- "mean"
         },
         var={ 
           thetaEst <- var(x)
           seThetaEst <- seVar(x)
           txtPar <- "variance"
         },
         dMean={
           thetaEst <- mean(x)-rho*mean(y)
           seThetaEst <- seDMean(x,y,rho=rho)
           txtPar <- paste("difference of ",ifelse(rho==1,"","(weighted) "),"means",sep="")
         },
         dVar={
           thetaEst <- var(x)-rho*var(y)
           seThetaEst <- seDVar(x,y,rho=rho)
           txtPar <- paste("difference of ",ifelse(rho==1,"","(weighted) "),"variances",sep="")
         },
         rMean={
           thetaEst <- mean(x)/mean(y)
           seThetaEst <- seRMean(x,y)
           txtPar <- "ratio of means"
         },
         rVar={
           thetaEst <- var(x)/var(y)
           seThetaEst <- seRVar(x,y)
           txtPar <- "ratio of variances"
         })
  estimate <- thetaEst
  names(estimate) <- txtPar
  statistic <- (thetaEst-reference)/seThetaEst
  names(statistic) <- "statistic"
  switch(alternative,
         less={ 
           p.value <-  pnorm(statistic)
           conf.int <- c(-Inf,thetaEst+seThetaEst*qnorm(conf.level))
         },
         greater={
           p.value <- 1-pnorm(statistic)
           conf.int <- c(thetaEst-seThetaEst*qnorm(conf.level),Inf)
         },
         two.sided={
           p.value <- 2*min(pnorm(statistic),1-pnorm(statistic))
           conf.int <- thetaEst+c(-1,1)*seThetaEst*qnorm(1- (1-conf.level)/2)
         })
  attr(conf.int,"conf.level") <- conf.level
  null.value <-  reference;names(null.value) <- txtPar		 
  if (parameter %in% c("mean","var")) tmp <- "One-sample" else tmp <- "Two-sample"
  method <- paste(tmp," asymptotic ",txtPar," test",sep="")
  
  res <-  list(statistic=statistic,p.value=p.value,conf.int=conf.int,
             estimate=estimate,null.value=null.value,alternative=alternative,method=method,data.name=data.name)
  class(res) <- "htest"
  res
  
}

asymp.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2)
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("asymp.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2)
        names(y$estimate) <- paste("mean in group", levels(g))
    y
}
