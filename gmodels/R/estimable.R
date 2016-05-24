## $Id: estimable.R 2060 2015-07-19 03:22:30Z warnes $
estimable <- function (obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
  {
    UseMethod("estimable")
  }

estimable.default <- function (obj, cm, beta0, conf.int=NULL,
                               show.beta0, joint.test=FALSE, ...)
{
  if (is.matrix(cm) || is.data.frame(cm))
    {
      cm <- t(apply(cm, 1, .to.est, obj=obj))
    }
  else if(is.list(cm))
    {
      cm <- matrix(.to.est(obj, cm), nrow=1)
    }
  else if(is.vector(cm))
    {
      cm <- matrix(.to.est(obj, cm), nrow=1)
    }
  else
    {
      stop("`cm' argument must be of type vector, list, or matrix.")
    }

  if(missing(show.beta0))
    {
      if(!missing(beta0))
        show.beta0=TRUE
      else
        show.beta0=FALSE
    }


  if (missing(beta0))
    {
      beta0 = rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm)))

    }


  if (joint.test==TRUE)
    {
      .wald(obj, cm, beta0)
    }
  else
    {
      if ("lme" %in% class(obj)) {
        stat.name <- "t.stat"
        cf <- summary(obj)$tTable
        rho <- summary(obj)$cor
        vcv <- rho * outer(cf[, 2], cf[, 2])
        tmp <- cm
        tmp[tmp==0] <- NA
        df.all <- t(abs(t(tmp) * obj$fixDF$X))
        df <- apply(df.all, 1, min, na.rm=TRUE)
        problem <- apply(df.all !=df, 1, any, na.rm=TRUE)
        if (any(problem))
          warning(paste("Degrees of freedom vary among parameters used to ",
                        "construct linear contrast(s): ",
                        paste((1:nrow(tmp))[problem], collapse=", "),
                        ". Using the smallest df among the set of parameters.",
                        sep=""))
      }
      else if ("lm" %in% class(obj))
        {
          stat.name <- "t.stat"
          cf <- summary.lm(obj)$coefficients
          vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
          df <- obj$df.residual
          if ("glm" %in% class(obj))
            {
              vcv <- summary(obj)$cov.scaled
              if(family(obj)[1] %in% c("poisson", "binomial"))
                {
                  stat.name <- "X2.stat"
                  df <- 1
                }
              else
                {
                  stat.name <- "t.stat"
                  df <- obj$df.residual
                }
            }
        }
      else if ("geese" %in% class(obj))
        {
          stat.name <- "X2.stat"
          cf <- summary(obj)$mean
          vcv <- obj$vbeta
          df <- 1
        }
      else if ("gee" %in% class(obj))
        {
          stat.name <- "X2.stat"
          cf <- summary(obj)$coef
          vcv <- obj$robust.variance
          df <- 1
        }
      else
        {
          stop("obj must be of class 'lm', 'glm', 'aov', 'lme', 'gee', 'geese' or 'nlme'")
        }
      if (is.null(cm))
        cm <- diag(dim(cf)[1])
      if (!dim(cm)[2]==dim(cf)[1])
        stop(paste("\n Dimension of ", deparse(substitute(cm)),
                   ": ", paste(dim(cm), collapse="x"),
                   ", not compatible with no of parameters in ",
                   deparse(substitute(obj)), ": ", dim(cf)[1], sep=""))
      ct <- cm %*% cf[, 1]
      ct.diff <- cm %*% cf[, 1] - beta0

      vc <- sqrt(diag(cm %*% vcv %*% t(cm)))
      if (is.null(rownames(cm)))
        rn <- paste("(", apply(cm, 1, paste, collapse=" "),
                    ")", sep="")
      else rn <- rownames(cm)
      switch(stat.name,
             t.stat={
               prob <- 2 * (1 - pt(abs(ct.diff/vc), df))
             },
             X2.stat={
               prob <- 1 - pchisq((ct.diff/vc)^2, df=1)
             })

      if (stat.name=="X2.stat")
        {
          retval <- cbind(hyp=beta0, est=ct, stderr=vc,
                          "X^2 value"=(ct.diff/vc)^2,
                          df=df, prob=1 - pchisq((ct.diff/vc)^2, df=1))
          dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
                                         "X^2 value", "DF", "Pr(>|X^2|)"))
        }
      else if (stat.name=="t.stat")
        {
          retval <- cbind(hyp=beta0, est=ct, stderr=vc, "t value"=ct.diff/vc,
                          df=df, prob=2 * (1 - pt(abs(ct.diff/vc), df)))
          dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
                                         "t value", "DF", "Pr(>|t|)"))
        }

      if (!is.null(conf.int))
        {
          if (conf.int <=0 || conf.int >=1)
            stop("conf.int should be between 0 and 1. Usual values are 0.95, 0.90")
          alpha <- 1 - conf.int
          switch(stat.name,
                 t.stat={
                   quant <- qt(1 - alpha/2, df)
                 },
                 X2.stat={
                   quant <- qt(1 - alpha/2, 100)
                 })
          nm <- c(colnames(retval), "Lower.CI", "Upper.CI")
          retval <- cbind(retval, lower=ct.diff - vc * quant, upper=ct.diff +
                          vc * quant)
          colnames(retval) <- nm
        }
      rownames(retval) <- make.unique(rownames(retval))
      retval <- as.data.frame(retval)
      if(!show.beta0) retval$beta0 <- NULL

      class(retval) <- c("estimable", class(retval))

      return(retval)
    }
}

.wald <- function (obj, cm,
                   beta0=rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm))))
{
  if (!is.matrix(cm) && !is.data.frame(cm))
    cm <- matrix(cm, nrow=1)
  df <- nrow(cm)
  if ("geese" %in% class(obj))
    {
      cf <- obj$beta
      vcv <- obj$vbeta
    }
  else if ("gee" %in% class(obj))
    {
      cf <- summary(obj)$coef
      vcv <- obj$robust.variance
    }
  else if ("lm" %in% class(obj))
    {
      cf <- summary.lm(obj)$coefficients[, 1]
      if ("glm" %in% class(obj))
        vcv <- summary(obj)$cov.scaled
      else
        vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
    }
  else if ("lme" %in% class(obj))
    {
      s.o <- summary(obj)
      cf <- s.o$tTable[,1]
      se <- s.o$tTable[, 2]
      rho <- s.o$cor
      vcv <- rho * outer(se, se)
    }
  else
    stop("obj must be of class 'lm', 'glm', 'aov', 'gee', 'geese', or 'lme'.")
  u <- (cm %*% cf)-beta0
  vcv.u <- cm %*% vcv %*% t(cm)
  W <- t(u) %*% solve(vcv.u) %*% u
  prob <- 1 - pchisq(W, df=df)
  retval <- as.data.frame(cbind(W, df, prob))
  names(retval) <- c("X2.stat", "DF", "Pr(>|X^2|)")
  print(as.data.frame(retval))
}

## estimable.mer <- function (obj, cm, beta0, conf.int=NULL, show.beta0,
##                            sim.mer=TRUE, n.sim=1000, ...)
## {
##   if (is.matrix(cm) || is.data.frame(cm))
##     {
##       cm <- t(apply(cm, 1, .to.est, obj=obj))
##     }
##   else if(is.list(cm))
##     {
##       cm <- matrix(.to.est(obj, cm), nrow=1)
##     }
##   else if(is.vector(cm))
##     {
##       cm <- matrix(.to.est(obj, cm), nrow=1)
##     }
##   else
##     {
##       stop("'cm' argument must be of type vector, list, or matrix.")
##     }

##   if(missing(show.beta0))
##     {
##       if(!missing(beta0))
##         show.beta0=TRUE
##       else
##         show.beta0=FALSE
##     }


##   if (missing(beta0))
##     {
##       beta0 = rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm)))

##     }

##   if ("mer" %in% class(obj)) {
##     if(sim.mer)
##       return(est.mer(obj=obj, cm=cm, beta0=beta0, conf.int=conf.int,
##                      show.beta0=show.beta0, n.sim=n.sim))

##     stat.name <- "mer"
##     cf <- as.matrix(fixef(obj))
##     vcv <- as.matrix(vcov(obj))
##     df <- NA
##   }
##   else {
##     stop("obj is not of class mer")
##   }

##   if (is.null(rownames(cm)))
##     rn <- paste("(", apply(cm, 1, paste, collapse=" "),
##                 ")", sep="")
##   else rn <- rownames(cm)

##   ct <- cm %*% cf[, 1]
##   ct.diff <- cm %*% cf[, 1] - beta0
##   vc <- sqrt(diag(cm %*% vcv %*% t(cm)))

##   retval <- cbind(hyp=beta0, est=ct, stderr=vc, "t value"=ct.diff/vc)
##   dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
##                                  "t value"))

##   rownames(retval) <- make.unique(rownames(retval))
##   retval <- as.data.frame(retval)
##   if(!show.beta0) retval$beta0 <- NULL

##   class(retval) <- c("estimable", class(retval))

##   return(retval)

## }
