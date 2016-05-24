irrelevant.args.warning <- function(object, args){
  if (any(!(names(object) %in% args))){
    irrelevant.args <- !(names(object) %in% args)
    be <- ifelse(sum(irrelevant.args) > 1, "are", "is")
    irrelevant.args <- paste(names(object)[irrelevant.args], collapse = ", ")
    warning(paste("arguments", irrelevant.args, be, "irrelevant and", be, "ignored", sep = " "))
  }
}

waldtest.mlogit <- function(object, ...){
  objects <- list(object, ...)
  margs <- c('nests', 'un.nest.el', 'unscaled', 'heterosc', 'rpar',
                   'R', 'correlation', 'halton', 'random.nb', 'panel')
  mlogit.args <- objects[names(objects) %in% margs]
  if (!is.null(names(objects))) objects <- objects[!(names(objects) %in% margs)]
  nmodels <- length(objects)
  specific.computation <- FALSE
  # if several models are provided, just use the default method
  if (nmodels > 1){
    return(waldtest.default(object, ...))
  }

  K <- length(colnames(model.matrix(object)))
  L <- length(object$freq)
  
  ###### guess the nature of the fitted model
  mixed.logit <- !is.null(object$call$rpar)
  heterosc.logit <- !is.null(object$call$heterosc) && object$call$heterosc
  nested.logit <- !is.null(object$call$nests)

  ###### Heteroscedastic logit model
  # the hypothesis is that J-1 parameters = 1
  if (heterosc.logit){
    su <- (K+1):(K+L-1)
    q <- rep(1, length(su))
    hyp <- "homoscedasticity"
  }
  
  ###### Nested logit Models
  if (nested.logit){
    J <- length(coef(object)) - K
    # First check whether the fitted model has a unique nest
    # elasticity or not
    if (is.null(object$call$un.nest.el)) un.nest.el <- FALSE
    else un.nest.el <- object$call$un.nest.el
    # If the fitted model has a unique nest elasticity, the only
    # relevant test is no nests : mlogit.args should be nests=NULL or
    # nothing. A warning is returned in case of supplementary arguments
    if (un.nest.el){
      if (!is.null(mlogit.args$nests)) stop("the nest argument should be NULL")
      irrelevant.args.warning(mlogit.args, "nests")
      su <- K + 1
      q <- 1
      hyp <- "no nests"
    }
    # If the nests elasticities are different, two possible tests :
    # 1. no nests (mlogit.args = (nests = NULL)) or nothing. stop if
    # !is.null(nests) and warning if other arguments than nests.
    # 2. unique nest elasticity (mlogit.args = (un.nest.el =
    # TRUE)). stop if un.nest.el = FALSE and warning if other arguments
    # are provided.
    if (!un.nest.el){
      if (!is.null(mlogit.args$nests)) stop("the nest argument should be NULL")
      if (!is.null(mlogit.args$un.nest.el) && mlogit.args$un.nest.el){
        irrelevant.args.warning(mlogit.args, "un.nest.el")
        su <- (K+1):length(coef(object))
        R <- matrix(0, nrow = length(coef(object)), ncol = length(su) - 1)
        for (i in 1:ncol(R)){
          R[K + 1, i] <- 1
          R[K + 1 + i, i] <- -1
        }
        Rb <- crossprod(R, coef(object))
        VRV <- t(R) %*% vcov(object) %*% R
        stat <- as.numeric(crossprod(Rb,solve(VRV, Rb)))
        df <- c(df = length(su) - 1)
        specific.computation <- TRUE
        hyp <- "unique nest elasticity"
      }
      else{
        if (length(mlogit.args) == 0 |
            ("nests" %in% names(mlogit.args) & is.null(mlogit.args$nests))){
          irrelevant.args.warning(mlogit.args, "nests")
          su <- (K+1):(K+J)
          q <- rep(1, length(su))
          hyp <- "no nests"
        }
        else{
          stop("irrelevant constrained model")
        }
      }
    }
  }

  ###### Mixed logit model
  if (mixed.logit){
    J <- length(object$rpar)
    # First check whether the random effects are correlated or not
    if (is.null(object$call$correlation)) correlation <- FALSE
    else correlation <- object$call$correlation
    # If the fitted model is uncorrelated, the only relevant test is
    # no random effects, mlogit.args = (rpar = NULL) ; stop if rpar is
    # not NULL and warning if supplementary arguments are provided
    if (!correlation){
      if (!is.null(mlogit.args$rpar)) stop("rpar should be NULL")
      irrelevant.args.warning(mlogit.args, "rpar")
      su <- K + (1:J)
      hyp <- "no random effects"
    }
    else{
      # if the fitted model is correlated, two possible tests :
      # 1. uncorrelated random effects : mlogit.args = (correlation =
      # FALSE), stop if (correlation = TRUE) and warning if
      # supplementary aguments are provided
      # 2. no random effects : mlogit.args = (rpar = NULL), stop if
      # rpar not NULL and a warning if supplementary arguments are
      # provided
      rd.el <- K+(1:(J*(J+1)/2))
      diag.el <- K + c(1, cumsum(J:2)+1)
      if (!is.null(mlogit.args$correlation) && mlogit.args$correlation)
        stop("irrelevant constrained model")
      if (!is.null(mlogit.args$correlation) && !mlogit.args$correlation){
        irrelevant.args.warning(mlogit.args, "correlation")
        su <- rd.el[!(rd.el %in% diag.el)]
        hyp <- "uncorrelated random effects"
      }
      else{
        if (!is.null(mlogit.args$rpar)) stop("rpar should be NULL")
        su <- rd.el
        hyp <- "no random effects"
      }
    }
    q <- rep(0, length(su))
  }
  
  if (!specific.computation){
    if (is.null(q)) wq <- coef(object)[su] else wq <- coef(object)[su] - q
    stat <- as.numeric(crossprod(wq,
                                 crossprod(solve(vcov(object)[su, su]),
                                           wq)))
    df <- c(df = length(su))
  }
  names(stat) <- 'chisq'
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  result <- list(statistic = stat,
                 parameter = df,
                 p.value = pval,
                 data.name = hyp,
                 method = "Wald test"
  #                 alternative = "unconstrainted model"
                 )
  class(result) <- 'htest'
  result
}

lrtest.mlogit <- function(object, ...){
  dots <- list(...)
  if (length(dots) == 0){
    model2 <- update(object, heterosc=FALSE, rpar = NULL,
                     start = NULL, nests = NULL,
                     gleontief = FALSE, method = 'nr')
    lrtest.default(object, model2)
  }
  else lrtest.default(object, ...)
}

scoretest <- function(object, ...){
  UseMethod("scoretest")
}

scoretest.mlogit <- function(object, ...){
  objects <- list(object, ...)
  margs <- c('nests', 'un.nest.el', 'unscaled', 'heterosc', 'rpar',
             'R', 'correlation', 'halton', 'random.nb', 'panel')
  mlogit.args <- objects[names(objects) %in% margs]
  if (!is.null(names(objects))) objects <- objects[!(names(objects) %in% margs)]
  nmodels <- length(objects)
  start.values <- c(coef(object))
  m <- list(nests = NULL, un.nest.el = FALSE, unscaled = FALSE, heterosc = FALSE,
            rpar = NULL, R = 40, correlation = FALSE, halton = NULL,
            random.nb = NULL, panel = FALSE)
  m[names(mlogit.args)] <- mlogit.args
  
  # if several models are provided, just use the default method
  if (nmodels > 1){
    return(scoretest.default(object, ...))
  }
  heterosc.logit <- (m$heterosc)
  nested.logit <- (!is.null(m$nests) || !is.null(object$nests))
  mixed.logit <- (!is.null(m$rpar) || m$correlation)
  if (heterosc.logit + nested.logit + mixed.logit == 0)
    stop("an unconstrained model should be described")
  if (heterosc.logit + nested.logit + mixed.logit > 1)
    stop("only one unconstrained model should be described")
  if (heterosc.logit){
    alt.hyp <- "heteroscedastic model"
    data.name <- "heterosc = TRUE"
  }
  if (nested.logit){
    init.nested.model <- !is.null(object$call$nests)
    if (init.nested.model){
      if (is.null(object$call$un.nest.el) || !object$call$un.nest.el){
        stop("irrelevant model for a score test")
      }
      J <- length(object$nests)
      start.values <- c(coef(object), rep(coef(object)[length(coef(object))], J - 1))
      data.name <- "un.nest.el = FALSE"
      alt.hyp <- "unique nest elasticity"
    }
    else{
      alt.hyp <- ifelse(m$un.nest.el, "nested model with a unique nest elasticity",
                        "nested model")
      nest.list <- c()
      for (i in 1:length(m$nests)){
        anest <- paste("c(\'",paste(m$nests[[i]],collapse="\',\'"),"\')", sep="")
        anest <- paste(names(m$nests)[i], " = ", anest, sep = "")
        nest.list <- c(nest.list, anest)
      }
      data.name = paste("nests = list(", paste(nest.list, collapse = ", "), ")", sep = "")
    }
  }
    
  if (mixed.logit){
    init.mixed.model <- !is.null(object$call$rpar)
    if (init.mixed.model){
      if (!is.null(object$call$correlation) && object$call$correlation) stop("not a relevant model for a score test")
      alt.hyp <- "uncorrelated random effects"
      data.name <- "correlation = TRUE"
    }
    else{
      if (m$correlation)
        alt.hyp <- "no correlated random effects"
      else alt.hyp <- "no uncorrelated random effects"
      data.name <- paste(names(m$rpar), paste("\'",as.character(m$rpar),"\'", sep = ""),
                         collapse = ",", sep = "=")
      data.name <- paste("rpar", "(", data.name, ")", sep = "")
    }
    if (init.mixed.model){
      J <- length(object$rpar)
      K <- ncol(model.matrix(object))
      sd <- coef(object)[-c(1:K)]
      rd.el <- K+(1:(J*(J+1)/2))
      diag.el <- K + c(1, cumsum(J:2)+1)
      start.values <- c(start.values[1:K], rep(0, length(rd.el)))
      start.values[diag.el] <- sd
    }
  }
  
  mc <- match.call()
  mc[[1]] <- as.name('update')
  mc[c('iterlim', 'method', 'start', 'print.level')] <- list(0, 'bfgs', start.values, 0)
  newmodel <- eval(mc, parent.frame())
  # gradient used to be a vector, now a matrix (the following ifelse should may be removed
  if (is.matrix(newmodel$gradient)) gradvect <- apply(newmodel$gradient, 2, sum) else gradvect <- newmodel$gradient
  stat <- - sum(gradvect * solve(newmodel$hessian, gradvect))
  names(stat) <- "chisq"
  df <- c(df = length(coef(newmodel)) - length(coef(object)))
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  result <- list(statistic = stat,
                 parameter = df,
                 p.value = pval,
                 data.name = data.name,
                 method = "score test",
                 alternative = alt.hyp
                 )
  class(result) <- 'htest'
  result
}

scoretest.default <- function(object, ...){
  new <- list(...)[[1]]
  cls <- class(object)[1]
  nmodels <- length(new)
  if (!inherits(new, 'formula') & !inherits(new, cls))
    stop("the updating argument doesn't have a correct class")
  if (inherits(new, cls)){
    ncoefs <- names(coef(new))
    new <- formula(formula(new))
  }
  else ncoefs <- names(coef(update(object, new, iterlim = 0)))
  start <- numeric(length = length(ncoefs))
  names(start) <- ncoefs
  supcoef <- ! ncoefs %in% names(coef(object))
  start[names(coef(object))] <- coef(object)
  newmodel <- update(object, new, start= start, iterlim = 0)
  data.name <- paste(deparse(formula(newmodel)))
  alt.hyp <- "unconstrained model"
  if (is.matrix(newmodel$gradient)) gradvect <- apply(newmodel$gradient, 2, sum) else gradvect <- newmodel$gradient
  stat <- - sum(gradvect * solve(newmodel$hessian, gradvect))
  names(stat) <- "chisq"
  df <- c(df = length(coef(newmodel)) - length(coef(object)))
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  result <- list(statistic = stat,
                 parameter = df,
                 p.value = pval,
                 data.name = data.name,
                 method = "score test",
                 alternative = alt.hyp
                 )
  class(result) <- 'htest'
  result
}
    
