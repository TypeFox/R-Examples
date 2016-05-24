setClass("FLXMRbetafix",
         representation(fixed_precision="ANY",
                        segment_precision="matrix",
                        design_precision="matrix"),
         contains = c("FLXMRbeta", "FLXMRfix"))

FLXMRbetafix <- function(formula = .~., precision = ~ 1,
                         fixed, fixed_precision, offset = NULL,
                         link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                         link.phi = "log", control = betareg.control())
{
  link <- match.arg(link)
  
  if (!is(precision, "formula"))
    stop("precision needs to be a formula")

  object <- FLXMRbeta(formula=formula, precision=precision, offset = offset,
                      link=link, link.phi=link.phi, control=control)
  if (!missing(fixed))
    object <- new("FLXMRbetafix", object, fixed=fixed, variance = FALSE)
  if (!missing(fixed_precision))
    object <- new("FLXMRbetafix", object, fixed_precision=fixed_precision, variance = TRUE)
  if (missing(fixed) && missing(fixed_precision)) stop("some parameters need to be fixed.")

  object@fit <- function(x, y, z, w, incidence, incidence_precision, segment, segment_precision, ...){
    fit <- betareg.fit(x, as.vector(y), z, weights=w, offset=offset, link = link, link.phi = link.phi,
                       control = control)
    k <- nrow(incidence)
    coefs <- coef(fit)
    names(coefs$mean) <- colnames(incidence)
    names(coefs$precision) <- colnames(incidence_precision)
    df <- rowSums(incidence/rep(colSums(incidence), each = nrow(incidence))) +
      rowSums(incidence_precision/rep(colSums(incidence_precision), each = nrow(incidence_precision))) 
    lapply(seq_len(k),
           function(K) with(list(coef = list(mean = coefs$mean[as.logical(incidence[K,])],
                                   precision = coefs$precision[as.logical(incidence_precision[K,])]),
                                 df = df[K],
                                 linkobjs = fit$link),
                            eval(object@defineComponent)))
  }
  object
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRbetafix"),
function(model, data, formula, lhs=TRUE, ...) {
  model_precision <- model
  model_precision@fullformula <- model_precision@formula <- . ~ .
  model_precision@fixed <- model@fixed_precision
  model <- getMethod("FLXgetModelmatrix", "FLXMRfix", where = asNamespace("flexmix"))(model, data, formula, lhs, ...)
  colnames(model@x) <- colnames(model@design)
  model_precision <- getMethod("FLXgetModelmatrix", "FLXMRfix", where = asNamespace("flexmix"))(model_precision, data, model@precision, lhs = FALSE, ...)
  model@z <- model_precision@x
  colnames(model@z) <- colnames(model_precision@design)
  model@terms_precision <- model_precision@terms
  model@contrasts_precision <- model_precision@contrasts
  model@xlevels_precision <- model_precision@xlevels
  model@design_precision <- model_precision@design
  model@segment_precision <- model_precision@segment
  model
})

setMethod("FLXmstep", signature(model = "FLXMRbetafix"), function(model, weights, ...)
{
  model@fit(model@x, model@y, model@z,
            as.vector(weights),
            model@design, model@design_precision,
            model@segment, model@segment_precision)
})

setMethod("predict", signature(object="FLXMRbetafix"),
function(object, newdata, components, ...)
{
  model <- FLXgetModelmatrix(object, newdata, object@fullformula, lhs=FALSE)
  k <- sum(object@nestedformula@k)
  N <- nrow(model@x)/k
  z <- list()
  for (m in seq_len(k)) {
    z[[m]] <- components[[m]]@predict(model@x[model@segment[,m], as.logical(model@design[m,]), drop=FALSE],
                                      model@z[model@segment_precision[,m], as.logical(model@design_precision[m,]), drop=FALSE], ...)
  }
  z
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRbetafix"), function(model, components, ...)
{
  sapply(seq_along(components), function(m)
         components[[m]]@logLik(model@x[model@segment[,m], as.logical(model@design[m,]), drop=FALSE],
                                model@y[model@segment[,m],,drop=FALSE],
                                model@z[model@segment_precision[,m], as.logical(model@design_precision[m,]), drop=FALSE]))
})

setMethod("FLXgetDesign", signature(object="FLXMRbetafix"),
function(object, ...) cbind(object@design, object@design_precision))

setMethod("FLXgetParameters", signature(object="FLXMRbetafix"),
function(object, components, ...) {
  parms <- lapply(components, function(x) unlist(slot(x, "parameters")))
  Design <- FLXgetDesign(object, components)
  parms_unique <- vector(length = ncol(Design))
  names(parms_unique) <- colnames(Design)
  for (k in seq_along(parms)) 
    parms_unique[as.logical(Design[k,])] <- parms[[k]]
  parms_unique
})

setMethod("FLXreplaceParameters", signature(object="FLXMRbetafix"),
function(object, components, parms) {
  Design <- list(mean = object@design, precision = object@design_precision)
  parms <- list(mean = parms[seq_len(ncol(Design$mean))],
                precision = parms[ncol(Design$mean) + seq_len(ncol(Design$precision))])
  lapply(seq_along(components), function(k) {
    Parameters <- list()
    for (i in seq_along(components[[k]]@parameters)) {
      parms_k <- parms[[i]][as.logical(Design[[i]][k,])]
      Parameters[[i]] <- parms_k[seq_along(components[[k]]@parameters[[i]])]
      attributes(Parameters[[i]]) <- attributes(components[[k]]@parameters[[i]])
      parms_k <- parms_k[-seq_along(components[[k]]@parameters[[i]])]
    }
    names(Parameters) <- names(components[[k]]@parameters)
    variables <- c("x", "y", "offset", "family")
    variables <- variables[variables %in% slotNames(object)]
    for (var in variables) 
      assign(var, slot(object, var))
    with(list(coef = Parameters, df = components[[k]]@df,
              linkobjs = get("linkobjs", environment(components[[k]]@predict))),
         eval(object@defineComponent))
  })
})

setMethod("FLXgradlogLikfun", signature(object="FLXMRbetafix"),
function(object, components, weights, ...) {
  Scores <- lapply(seq_along(components), function(k) {
    linkobjs <- get("linkobjs", environment(components[[k]]@predict))
    pp <- components[[k]]@predict(object@x[object@segment[,k], as.logical(object@design[k,]), drop=FALSE],
                                  object@z[object@segment_precision[,k], as.logical(object@design_precision[k,]), drop=FALSE], ...)
    ystar <- qlogis(as.vector(object@y[object@segment[,k]]))
    mu <- pp$mean
    phi <- pp$precision
    mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
    res <- phi * (ystar - mustar) * linkobjs$mean$mu.eta(linkobjs$mean$linkfun(mu))
    Scores_x <- weights[,k] * res * object@x[object@segment[,k], as.logical(object@design[k,]), drop=FALSE]
    res <- (mu * (ystar - mustar) + log(1 - as.vector(object@y[object@segment[,k]])) -
            digamma((1-mu) * phi) + digamma(phi)) * linkobjs$precision$mu.eta(linkobjs$precision$linkfun(phi))
    Scores_z <- weights[,k] * res * object@z[object@segment_precision[,k], as.logical(object@design_precision[k,]), drop=FALSE]
    new_x <- (colSums(object@design[seq_len(k),,drop=FALSE]) == 1)[as.logical(object@design[k,])]
    new_z <- (colSums(object@design_precision[seq_len(k),,drop=FALSE]) == 1)[as.logical(object@design_precision[k,])]
    list(mean = Scores_x[,new_x,drop=FALSE], precision = Scores_z[,new_z,drop=FALSE])
  })
  lapply(Scores, function(x) do.call("cbind", x))
})                       

