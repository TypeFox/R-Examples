effectsplot_psychomix <- function(object,
  ask = FALSE, confint = FALSE, style = c("lines", "stacked"), colors = NULL, ...)
{
  if(!(inherits(object, "efflist") | inherits(object, "effpoly") | inherits(object, "eff"))) object <- effects::allEffects(object)
  if(inherits(object, "efflist")) {
    stopifnot(inherits(object[[1]], "effpoly") | inherits(object[[1]], "eff"))
    k <- length(object[[1]]$y.levels)
    for(i in seq_along(object)) object[[i]]$response <- "Component"
  } else {
    k <- length(object$y.levels)
    object$response <- "Component"
  }
  style <- match.arg(style)
  if(is.null(colors)) colors <- if(style == "stacked") gray.colors(k) else qualitative_hcl(k)
  return(plot(object, ask = ask, confint = confint, style = style, colors = colors, ...))
}

refit_concomitant_psychomix <- function(object, binary = c("glm", "multinom"), ...)
{
  ## classes/weights
  p <- posterior(object)
  k <- ncol(p)
  .weights <- as.vector(p)

  ## check for concomitants
  f <- as.formula(object@call$formula)
  tm <- f <- delete.response(terms(f)) ## FIXME: but need to handle extreme scorers
  if(length(tm) < 1L) stop("no concomitant variables")
  attributes(f) <- NULL
  f <- as.formula(f)
  
  ## get original model frame
  mf <- object@call
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- f
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## expand and add auxiliary response/weights
  n <- nrow(mf)
  mf <- mf[rep(1:n, k), , drop = FALSE]
  mf$.Component <- factor(rep(1:k, each = n))
  mf$.weights <- .weights
  
  ## refit concomitant model
  f <- update(f, .Component ~ .)
  
  ## effects package cannot treat multinom object with binary response
  if(k == 2) {
    if(match.arg(binary) == "glm")
      return(suppressWarnings(glm(f, data = mf, weights = .weights, family = binomial)))
  }
  nnet::multinom(f, data = mf, weights = .weights, trace = FALSE)
}

setMethod("effectsplot", "btmix", effectsplot_psychomix)
setMethod("effectsplot", "raschmix", effectsplot_psychomix)
effectsplot.efflist <- effectsplot_psychomix
effectsplot.effpoly <- effectsplot_psychomix
effectsplot.eff     <- effectsplot_psychomix

effect.raschmix <- effect.btmix <- function(term, mod, ...) effects::effect(term, refit_concomitant_psychomix(mod), ...)
allEffects.raschmix <- allEffects.btmix <- function(object, ...) effects::allEffects(refit_concomitant_psychomix(object), ...)

