#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: flexmix.R 4978 2014-02-13 15:45:15Z gruen $
#

log_row_sums <- function(m) {
  M <- m[cbind(seq_len(nrow(m)), max.col(m))]
  M + log(rowSums(exp(m - M)))
}

## The following two methods only fill in and rearrange the model argument
setMethod("flexmix",
          signature(formula = "formula", model="missing"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
    mycall = match.call()
    z <- flexmix(formula=formula, data=data, k=k, cluster=cluster,
                 model=list(FLXMRglm()), concomitant=concomitant,
                 control=control, weights = weights)
    z@call <- mycall
    z
})

setMethod("flexmix",
          signature(formula = "formula", model="FLXM"),
function(formula, data=list(), k=NULL, cluster=NULL, 
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
  mycall = match.call()
  z <- flexmix(formula=formula, data=data, k=k, cluster=cluster, 
               model=list(model), concomitant=concomitant,
               control=control, weights=weights)
  z@call <- mycall
  z
})


## This is the real thing
setMethod("flexmix",
          signature(formula = "formula", model="list"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
    mycall = match.call()
    control = as(control, "FLXcontrol")
    if (!is(concomitant, "FLXP")) concomitant <- FLXPconstant()
    
    groups <- .FLXgetGrouping(formula, data)
    model <- lapply(model, FLXcheckComponent, k, cluster)
    k <- unique(unlist(sapply(model, FLXgetK, k)))
    if (length(k) > 1) stop("number of clusters not specified correctly")

    model <- lapply(model, FLXgetModelmatrix, data, formula)
    
    groups$groupfirst <-
        if (length(groups$group)) groupFirst(groups$group)
        else rep(TRUE, FLXgetObs(model[[1]]))
    
    if (is(weights, "formula")) {
      weights <- model.frame(weights, data = data, na.action = NULL)[,1]
    }
    ## check if the weights are integer
    ## if non-integer weights are wanted modifications e.g.
    ## for classify != weighted and
    ## plot,flexmix,missing-method are needed
    if (!is.null(weights) & !identical(weights, as.integer(weights)))
      stop("only integer weights allowed")
    ## if weights and grouping is specified the weights within each
    ## group need to be the same
    if (!is.null(weights) & length(groups$group)>0) {
      unequal <- tapply(weights, groups$group, function(x) length(unique(x)) > 1)
      if (any(unequal)) stop("identical weights within groups needed")
    }
    
    postunscaled <- initPosteriors(k, cluster, FLXgetObs(model[[1]]), groups)

    if (ncol(postunscaled) == 1L)
      concomitant <- FLXPconstant()
  
    concomitant <- FLXgetModelmatrix(concomitant, data = data,
                                     groups = groups)
    

    z <- FLXfit(model=model, concomitant=concomitant, control=control,
                postunscaled=postunscaled, groups=groups, weights = weights)
    
    z@formula = formula
    z@call = mycall
    z@k0 = as.integer(k)
    z
})

###**********************************************************

setMethod("FLXgetK", signature(model = "FLXM"), function(model, k, ...) k)
setMethod("FLXgetObs", signature(model = "FLXM"), function(model) nrow(model@x))
setMethod("FLXcheckComponent", signature(model = "FLXM"), function(model, ...) model)
setMethod("FLXremoveComponent", signature(model = "FLXM"), function(model, ...) model)

setMethod("FLXmstep", signature(model = "FLXM"), function(model, weights, components, ...) {
  if ("component" %in% names(formals(model@fit)))
    sapply(seq_len(ncol(weights)), function(k) model@fit(model@x, model@y, weights[,k], component = components[[k]]@parameters))
  else
    sapply(seq_len(ncol(weights)), function(k) model@fit(model@x, model@y, weights[,k]))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXM"), function(model, components, ...) {
  matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
})

###**********************************************************
setMethod("FLXfit", signature(model="list"),
function(model, concomitant, control, postunscaled=NULL, groups, weights)
{
  ### initialize
  k <- ncol(postunscaled)
  N <- nrow(postunscaled)
  control <- allweighted(model, control, weights)
  if(control@verbose>0)
    cat("Classification:", control@classify, "\n")
  if (control@classify=="random") iter.rm <- 0
  group <- groups$group
  groupfirst <- groups$groupfirst
  if(length(group)>0) postunscaled <- groupPosteriors(postunscaled, group)

  logpostunscaled <- log(postunscaled)
  postscaled <- exp(logpostunscaled - log_row_sums(logpostunscaled))
  
  llh <- -Inf
  if (control@classify=="random") llh.max <- -Inf
  converged <- FALSE
  components <- rep(list(rep(list(new("FLXcomponent")), k)), length(model))
  ### EM
  for(iter in seq_len(control@iter.max)) {
      ### M-Step
      postscaled = .FLXgetOK(postscaled, control, weights)
      prior <- if (is.null(weights))
        ungroupPriors(concomitant@fit(concomitant@x, postscaled[groupfirst,,drop=FALSE]),
                      group, groupfirst)
      else ungroupPriors(concomitant@fit(concomitant@x, (postscaled/weights)[groupfirst & weights > 0,,drop=FALSE], weights[groupfirst & weights > 0]),
                         group, groupfirst)
      # Check min.prior
      nok <- if (nrow(prior) == 1) which(prior < control@minprior) else {
               if (is.null(weights)) which(colMeans(prior[groupfirst,,drop=FALSE]) < control@minprior)
               else which(colSums(prior[groupfirst,] * weights[groupfirst])/sum(weights[groupfirst]) < control@minprior)
             }
      if(length(nok)) {
        if(control@verbose>0)
          cat("*** Removing", length(nok), "component(s) ***\n")
        prior <- prior[,-nok,drop=FALSE]
        prior <- prior/rowSums(prior)
        postscaled <- postscaled[,-nok,drop=FALSE]
        postscaled[rowSums(postscaled) == 0,] <- if (nrow(prior) > 1) prior[rowSums(postscaled) == 0,]
                                                 else prior[rep(1, sum(rowSums(postscaled) == 0)),]
        postscaled <- postscaled/rowSums(postscaled)
        if (!is.null(weights)) postscaled <- postscaled * weights
        k <- ncol(prior)
        if (k == 0) stop("all components removed")
        if (control@classify=="random") {
          llh.max <- -Inf
          iter.rm <- iter
        }
        model <- lapply(model, FLXremoveComponent, nok)
      }
      components <- lapply(seq_along(model), function(i) FLXmstep(model[[i]], postscaled, components[[i]]))
      postunscaled <- matrix(0, nrow = N, ncol = k)
      for (n in seq_along(model))
        postunscaled <- postunscaled + FLXdeterminePostunscaled(model[[n]], components[[n]])
      if(length(group)>0)
        postunscaled <- groupPosteriors(postunscaled, group)
      ### E-Step
      ## Code changed thanks to Nicolas Picard
      ## to avoid problems with small likelihoods
      postunscaled <- if (nrow(prior) > 1) postunscaled + log(prior)
                         else sweep(postunscaled, 2, log(prior), "+")
      logpostunscaled <- postunscaled
      postunscaled <- exp(postunscaled)
      postscaled <- exp(logpostunscaled - log_row_sums(logpostunscaled))
      ##<FIXME>: wenn eine beobachtung in allen Komonenten extrem
      ## kleine postunscaled-werte hat, ist exp(-postunscaled)
      ## numerisch Null, und damit postscaled NaN
      ## log(rowSums(postunscaled)) ist -Inf
      ##</FIXME>
      if (any(is.nan(postscaled))) {
        index <- which(as.logical(rowSums(is.nan(postscaled))))
        postscaled[index,] <- if(nrow(prior)==1) rep(prior, each = length(index)) else prior[index,]
        postunscaled[index,] <- .Machine$double.xmin
      }
      ### check convergence
      llh.old <- llh
      llh <- if (is.null(weights)) sum(log_row_sums(logpostunscaled[groupfirst,,drop=FALSE]))
             else sum(log_row_sums(logpostunscaled[groupfirst,,drop=FALSE])*weights[groupfirst])
      if(is.na(llh) | is.infinite(llh))
        stop(paste(formatC(iter, width=4),
                   "Log-likelihood:", llh))
      if (abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
        if(control@verbose>0){
          printIter(iter, llh)
          cat("converged\n")
        }
        converged <- TRUE
        break
      }
      if (control@classify=="random") {
        if (llh.max < llh) {
          components.max <- components
          prior.max <- prior
          postscaled.max <- postscaled
          postunscaled.max <- postunscaled
          llh.max <- llh
        }
      }
      if(control@verbose && (iter%%control@verbose==0))
        printIter(iter, llh)
    }
  ### Construct return object
  if (control@classify=="random") {
    components <- components.max
    prior <- prior.max
    postscaled <- postscaled.max
    postunscaled <- postunscaled.max
    llh <- llh.max
    iter <- control@iter.max - iter.rm
  }

  components <- lapply(seq_len(k), function(i) lapply(components, function(x) x[[i]]))
  names(components) <- paste("Comp", seq_len(k), sep=".") 
  cluster <- max.col(postscaled)
  size <-  if (is.null(weights)) tabulate(cluster, nbins=k) else tabulate(rep(cluster, weights), nbins=k)
  names(size) <- seq_len(k)
  concomitant <- FLXfillConcomitant(concomitant, postscaled[groupfirst,,drop=FALSE], weights[groupfirst])
  df <- concomitant@df(concomitant@x, k) + sum(sapply(components, sapply, slot, "df"))
  control@nrep <- 1
  prior <- if (is.null(weights)) colMeans(postscaled[groupfirst,,drop=FALSE])
           else colSums(postscaled[groupfirst,,drop=FALSE] * weights[groupfirst])/sum(weights[groupfirst])

  retval <- new("flexmix", model=model, prior=prior,
                posterior=list(scaled=postscaled,
                  unscaled=postunscaled),
                weights = weights,
                iter=iter, cluster=cluster, size = size,
                logLik=llh, components=components,
                concomitant=concomitant,
                control=control, df=df, group=group, k=as(k, "integer"),
                converged=converged)
  retval
})

###**********************************************************
.FLXgetOK = function(p, control, weights){
    n = ncol(p)
    N = seq_len(n)
    if (is.null(weights)) {
      if (control@classify == "weighted")
        return(p)
      else {
        z = matrix(FALSE, nrow = nrow(p), ncol = n)
        if(control@classify %in% c("CEM", "hard")) 
            m = max.col(p)
        else if(control@classify %in% c("SEM", "random")) 
            m = apply(p, 1, function(x) sample(N, size = 1, prob = x))
        else stop("Unknown classification method")
        z[cbind(seq_len(nrow(p)), m)] = TRUE
      }
    }else {
      if(control@classify=="weighted")
        z <- p * weights
      else{
        z = matrix(FALSE,  nrow=nrow(p), ncol=n)       
        if(control@classify %in% c("CEM", "hard")) {
          m = max.col(p)
          z[cbind(seq_len(nrow(p)), m)] = TRUE
          z <- z * weights
        }
        else if(control@classify %in% c("SEM", "random")) 
          z = t(sapply(seq_len(nrow(p)), function(i) table(factor(sample(N, size=weights[i], prob=p[i,], replace=TRUE), N))))
        else stop("Unknown classification method")
      }
    }
    z
}    

###**********************************************************

RemoveGrouping <- function(formula) {
  lf <- length(formula)
  formula1 <- formula
  if(length(formula[[lf]])>1) {
    if (deparse(formula[[lf]][[1]]) == "|"){
      formula1[[lf]] <- formula[[lf]][[2]]
    }
    else if (deparse(formula[[lf]][[1]]) == "("){
      form <- formula[[lf]][[2]]
      if (length(form) == 3 && form[[1]] == "|")
        formula1[[lf]] <- form[[2]]
    }
  }
  formula1
}

.FLXgetGroupingVar <- function(x)
{
  lf <- length(x)
  while (lf > 1) {
    x <- x[[lf]]
    lf <- length(x)
  }
  x
}

.FLXgetGrouping <- function(formula, data)
{
  group <- factor(integer(0))
  formula1 <- RemoveGrouping(formula)
  if (!identical(formula1, formula))
    group <- factor(eval(.FLXgetGroupingVar(formula), data))
  return(list(group=group, formula=formula1))
}

setMethod("FLXgetModelmatrix", signature(model="FLXM"),
function(model, data, formula, lhs=TRUE, ...)
{
  formula <- RemoveGrouping(formula)
  if (length(grep("\\|", deparse(model@formula)))) stop("no grouping variable allowed in the model")
  if(is.null(model@formula))
    model@formula = formula
  
  ## model@fullformula = update.formula(formula, model@formula)
  ## <FIXME>: ist das der richtige weg, wenn ein punkt in beiden
  ## formeln ist?
  model@fullformula = update(terms(formula, data=data), model@formula)
  ## </FIXME>
  
  if (lhs) {
    mf <- if (is.null(model@terms)) model.frame(model@fullformula, data=data, na.action = NULL) else model.frame(model@terms, data=data, na.action = NULL, xlev = model@xlevels)
    model@terms <- attr(mf, "terms")
    response <- as.matrix(model.response(mf))
    model@y <- model@preproc.y(response)
  }
  else {
    mt1 <- if (is.null(model@terms)) terms(model@fullformula, data=data) else model@terms
    mf <- model.frame(delete.response(mt1), data=data, na.action = NULL, xlev = model@xlevels)
    model@terms<- attr(mf, "terms")
    ## <FIXME>: warum war das da???
    ## attr(mt, "intercept") <- attr(mt1, "intercept")
    ## </FIXME>
  }
  X <- model.matrix(model@terms, data=mf)
  model@contrasts <- attr(X, "contrasts")
  model@x <- model@preproc.x(X)
  model@xlevels <- .getXlevels(model@terms, mf)
  model
})

## groupfirst: for grouped observation we need to be able to use
## the posterior of each group, but for computational simplicity
## post(un)scaled has N rows (with mutiple identical rows for each
## group). postscaled[groupfirst,] extracts posteriors of each
## group ordered as the appear in the data set.
groupFirst <- function(x) !duplicated(x)

## if we have a group variable, set the posterior to the product
## of all density values for that group (=sum in logarithm)
groupPosteriors <- function(x, group)
{    
    for(g in levels(group)){
        gok <- group==g
        if(any(gok)){
            x[gok,] <- matrix(colSums(x[gok,,drop=FALSE]),
                              nrow=sum(gok), ncol=ncol(x), byrow=TRUE)
        }
    }
    x
}

ungroupPriors <- function(x, group, groupfirst) {
  if (!length(group)) group <- seq_along(groupfirst)
  if (nrow(x) >= length(group[groupfirst])) {
    x <- x[order(as.integer(group[groupfirst])),,drop=FALSE]
    x <- x[as.integer(group),,drop=FALSE]
  }
  x
}

allweighted <- function(model, control, weights) {
  allweighted <- all(sapply(model, function(x) x@weighted))
  if(allweighted){
    if(control@classify=="auto")
      control@classify <- "weighted"
  }
  else{
    if(control@classify=="auto") 
      control@classify <- "hard"
    else if (control@classify=="weighted") {
      warning("only hard classification supported for the modeldrivers")
      control@classify <- "hard"
    }
      
    if(!is.null(weights))
      stop("it is not possible to specify weights for models without weighted ML estimation")
  }
  control
}

initPosteriors <- function(k, cluster, N, groups) {
  if(is(cluster, "matrix")){
    postunscaled <- cluster
    if (!is.null(k)) if (k != ncol(postunscaled)) stop("specified k does not match the number of columns of cluster")
  }
  else{
    if(is.null(cluster)){
      if(is.null(k))
        stop("either k or cluster must be specified")
      else
        cluster <- ungroupPriors(as.matrix(sample(seq_len(k), size = sum(groups$groupfirst), replace=TRUE)),
                                           groups$group, groups$groupfirst)
    }
    else{
      cluster <- as(cluster, "integer")
      if (!is.null(k)) if (k != max(cluster)) stop("specified k does not match the values in cluster")
      k <- max(cluster)
    }
    postunscaled <- matrix(0.1, nrow=N, ncol=k)
    for(K in seq_len(k)){
      postunscaled[cluster==K, K] <- 0.9
    }
  }
  postunscaled
}

###**********************************************************

setMethod("predict", signature(object="FLXdist"),
function(object, newdata=list(), aggregate=FALSE, ...){
    if (missing(newdata)) return(fitted(object, aggregate=aggregate, drop=FALSE))
    x = list()
    for(m in seq_along(object@model)) {
      comp <- lapply(object@components, "[[", m)
      x[[m]] <- predict(object@model[[m]], newdata, comp, ...)
    }
    if (aggregate) {
      prior_weights <- prior(object, newdata)
      z <- lapply(x, function(z) matrix(rowSums(do.call("cbind", z) * prior_weights), nrow = nrow(z[[1]])))
    }
    else {
      z <- list()
      for (k in seq_len(object@k)) {
        z[[k]] <- do.call("cbind", lapply(x, "[[", k))
      }
      names(z) <- paste("Comp", seq_len(object@k), sep=".")
    }
    z
})

###**********************************************************

setMethod("parameters", signature(object="FLXdist"),
function(object, component=NULL, model=NULL, which = c("model", "concomitant"),
         simplify=TRUE, drop=TRUE)
{
    which <- match.arg(which)
    if (is.null(component)) component <- seq_len(object@k)
    if (is.null(model)) model <- seq_along(object@model)

    if (which == "model") {
      if (simplify) {
        parameters <- sapply(model, function(m) sapply(object@components[component], function(x) unlist(x[[m]]@parameters), simplify=TRUE),
                             simplify = FALSE)
      }
      else {
        parameters <- sapply(model, function(m) sapply(object@components[component], function(x) x[[m]]@parameters, simplify=FALSE),
                             simplify = FALSE)
      }
      if (drop) {
        if (length(component) == 1 && !simplify) parameters <- lapply(parameters, "[[", 1)
        if (length(model) == 1) parameters <- parameters[[1]]
      }
    } else {
      parameters <- object@concomitant@coef[,component,drop=FALSE]
    }
    parameters
})

setMethod("prior", signature(object="FLXdist"),
function(object, newdata, ...) {
  if (missing(newdata))
    prior <- object@prior
  else {
    groups <- .FLXgetGrouping(object@formula, newdata)
    nobs <- if (is(newdata, "data.frame")) nrow(newdata)
            else min(sapply(newdata, function(x) {
              if (is(x, "matrix")) nrow(x) else length(x)
            }))
    group <- if (length(groups$group)) groups$group else factor(seq_len(nobs))
    object@concomitant <- FLXgetModelmatrix(object@concomitant, data = newdata,
                                            groups = list(group=group,
                                              groupfirst = groupFirst(group)))
    prior <- determinePrior(object@prior, object@concomitant, group)[as.integer(group),]
  }
  prior
})


setMethod("posterior", signature(object="flexmix", newdata="missing"),
function(object, newdata, unscaled = FALSE, ...)
{
  if (unscaled) return(object@posterior$unscaled)
  else return(object@posterior$scaled)
})

setMethod("posterior", signature(object="FLXdist", newdata="listOrdata.frame"),
          function(object, newdata, unscaled=FALSE,...) {
            comp <- lapply(object@components, "[[", 1)
            postunscaled <- posterior(object@model[[1]], newdata, comp, ...)
            for (m in seq_along(object@model)[-1]) {
              comp <- lapply(object@components, "[[", m)
              postunscaled <- postunscaled + posterior(object@model[[m]], newdata, comp, 
                                                       ...)
            }
            groups <- .FLXgetGrouping(object@formula, newdata)
            prior <- prior(object, newdata = newdata)
            if(length(groups$group)>0)
              postunscaled <- groupPosteriors(postunscaled, groups$group)
            postunscaled <- postunscaled + log(prior)
            if (unscaled) return(exp(postunscaled))
            else return(exp(postunscaled - log_row_sums(postunscaled)))
})            

setMethod("posterior", signature(object="FLXM", newdata="listOrdata.frame"),
function(object, newdata, components, ...) {
  object <- FLXgetModelmatrix(object, newdata, object@fullformula, lhs = TRUE)
  FLXdeterminePostunscaled(object, components, ...)
})
    
setMethod("clusters", signature(object="flexmix", newdata="missing"),
function(object, newdata, ...)
{
    object@cluster
})
    
setMethod("clusters", signature(object="FLXdist", newdata="ANY"),
function(object, newdata, ...)
{
    max.col(posterior(object, newdata, ...))
})

###**********************************************************

setMethod("summary", "flexmix",
function(object, eps=1e-4, ...){    
    z <- new("summary.flexmix",
             call = object@call,
             AIC = AIC(object),
             BIC = BIC(object),
             logLik = logLik(object))

    TAB <- data.frame(prior=object@prior,
                      size=object@size)
    rownames(TAB) <- paste("Comp.", seq_len(nrow(TAB)), sep="")
    TAB[["post>0"]] <- if (is.null(object@weights)) colSums(object@posterior$scaled > eps)
                       else colSums((object@posterior$scaled > eps) * object@weights)
    TAB[["ratio"]] <- TAB[["size"]]/TAB[["post>0"]]
    
    z@comptab = TAB
    z
    
})

###**********************************************************

