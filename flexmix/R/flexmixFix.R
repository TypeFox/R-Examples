#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: flexmixFix.R 4978 2014-02-13 15:45:15Z gruen $
#

setMethod("FLXcheckComponent", signature(model = "FLXMRfix"), function(model, k, cluster, ...) {
  if (sum(model@nestedformula@k)) {
    if (!is.null(k)) {
      if (k != sum(model@nestedformula@k)) stop("specified k does not match the nestedformula in the model")
    }
    else k <- sum(model@nestedformula@k)
  }
  else {
    if (is(cluster, "matrix")) {
      if (is.null(k)) k <- ncol(cluster)
    }
    else if (!is.null(cluster)) {
      if (is.null(k)) {
        cluster <- as(cluster, "integer")
        k <- max(cluster)
      }
    }
    if (is.null(k)) stop("either k, cluster or the nestedformula of the model must be specified")
    else model@nestedformula <- as(k, "FLXnested")
  }
  if (length(model@variance) > 1) {
    if (sum(model@variance) != k) stop("specified k does not match the specified varFix argument in the model")
  }
  else if (model@variance) model@variance <- k
  else model@variance <- rep(1, k)

  model
})

setMethod("FLXgetObs", signature(model = "FLXMRfix"), function(model) nrow(model@y)/sum(model@nestedformula@k))
setMethod("FLXgetK", signature(model = "FLXMRfix"), function(model, ...) sum(model@nestedformula@k))

setMethod("FLXremoveComponent", signature(model = "FLXMRfix"), function(model, nok, ...)
{
  if (!length(nok)) return(model)
  K <- model@nestedformula
  wnok <- sapply(nok, function(i) which(apply(rbind(i > c(0, cumsum(K@k[-length(K@k)])),
                                                    i <= c(cumsum(K@k))), 2, all)))
  for (w in wnok) {
    K@k[w] <- K@k[w] - 1
    if (K@k[w] == 0) {
      K@k <- K@k[-w]
      K@formula <- K@formula[-w]
    }
    k <- sum(K@k)
    model@nestedformula <- K
  }
  varnok <- sapply(nok, function(i) which(apply(rbind(i > c(0, cumsum(model@variance[-length(model@variance)])),
                                                      i <= c(cumsum(model@variance))), 2, all)))
  for (w in varnok) {
    model@variance[w] <- model@variance[w] - 1
    if (model@variance[w] == 0) {
      model@variance <- model@variance[-w]
    }
  }
  
  rok <- which(!apply(model@segment[,nok,drop=FALSE], 1, function(x) any(x)))
  model@x <- model@x[rok, which(colSums(model@design[-nok,,drop=FALSE]) > 0), drop=FALSE]
  model@y <- model@y[rok,, drop=FALSE]
  model@design <- model@design[-nok,,drop=FALSE]
  cok <- colSums(model@design) > 0
  model@design <- model@design[,cok,drop=FALSE]
  model@segment <- model@segment[rok,-nok, drop=FALSE]
  model
})

###**********************************************************
setMethod("FLXmstep", signature(model = "FLXMRfix"), function(model, weights, ...)
{
  model@fit(model@x, model@y,
            as.vector(weights),
            model@design, model@variance)
})

###**********************************************************
setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRfix"), function(model, components, ...)
{
  sapply(seq_along(components), function(m)
         components[[m]]@logLik(model@x[model@segment[,m], as.logical(model@design[m,]), drop=FALSE],
                                model@y[model@segment[,m],,drop=FALSE]))
})

###**********************************************************

modelMatrix <- function(random, fixed, nested, data=list(), lhs, xlevels = NULL)
{
  if (!lhs) 
    random <- random[-2]
  mf.random <- model.frame(random, data=data, na.action = NULL)
  response <- if (lhs) as.matrix(model.response(mf.random)) else NULL
  xlev <- xlevels[names(.getXlevels(terms(mf.random), mf.random))]
  mm.random <- if (is.null(xlev))  model.matrix(terms(mf.random), data=mf.random)
               else model.matrix(terms(mf.random), data=data, xlev = xlev)
  xlevels.random <- .getXlevels(terms(mf.random), mf.random)
  randomfixed <- if(identical(paste(deparse(fixed), collapse = ""), "~0")) random
                 else update(random, paste("~.+", paste(deparse(fixed[[length(fixed)]]), collapse = "")))
  mf.randomfixed <- model.frame(randomfixed, data=data)
  mm.randomfixed <- model.matrix(terms(mf.randomfixed), data=mf.randomfixed, xlev = xlevels[names(.getXlevels(terms(mf.randomfixed), mf))])
  mm.fixed <- mm.randomfixed[,!colnames(mm.randomfixed) %in% colnames(mm.random), drop=FALSE]
  xlevels.fixed <- .getXlevels(terms(mf.randomfixed), mf.randomfixed)
  all <- mm.all <- mm.nested <- xlevels.nested <- list()
  for (l in seq_along(nested)) {
    all[[l]] <- if (identical(paste(deparse(nested[[l]]), collapse = ""), "~0")) randomfixed
                else update(randomfixed, paste("~.+", paste(deparse(nested[[l]][[length(nested[[l]])]]), collapse = "")))
    mf <- model.frame(all[[l]], data=data)
    mm.all[[l]] <- model.matrix(terms(mf), data=mf, xlev = xlevels[names(.getXlevels(terms(mf), mf))])
    mm.nested[[l]] <- mm.all[[l]][,!colnames(mm.all[[l]]) %in% colnames(mm.randomfixed),drop=FALSE]
    xlevels.nested[[l]] <- .getXlevels(terms(mf), mf)
  }
  return(list(random=mm.random, fixed=mm.fixed, nested=mm.nested, response=response, xlevels=c(xlevels.random, xlevels.fixed, unlist(xlevels.nested))))
}

###**********************************************************

modelDesign <- function(mm.all, k) {
  design <- matrix(1, nrow=sum(k@k), ncol=ncol(mm.all$fixed))
  col.names <- colnames(mm.all$fixed)
  nested <- matrix(0, nrow=sum(k@k), ncol=sum(sapply(mm.all$nested, ncol)))
  cumK <- c(0, cumsum(k@k))
  i <- 0
  for (l in seq_along(mm.all$nested)) {
    if (ncol(mm.all$nested[[l]])) {
      nested[(cumK[l] + 1):cumK[l+1], i+seq_len(ncol(mm.all$nested[[l]]))] <- 1
      i <- i+ncol(mm.all$nested[[l]])
      col.names <- c(col.names, colnames(mm.all$nested[[l]]))
    }
  }
  design <- cbind(design, nested)
  if (ncol(mm.all$random)) design <- cbind(design,
                                           kronecker(diag(sum(k@k)), matrix(1, ncol=ncol(mm.all$random))))
  colnames(design) <- c(col.names, rep(colnames(mm.all$random), sum(k@k)))
  design
}

###**********************************************************



