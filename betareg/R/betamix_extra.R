setClass("FLXMRbeta_extra",
         representation(extra_components="ANY"),
         contains = "FLXMRbeta")

extraComponent <- function(type = c("uniform", "betareg"), coef, delta, link = "logit", link.phi = "log") {
  type <- match.arg(type)
  z <- list(coef = coef)
  class(z) <- "extraComponent"
  if (type == "uniform") {
    if (missing(delta)) stop("specify a 'delta' for the half-length of the uniform distribution")
    else z$delta <- delta
  }
  if (type == "betareg") {
    if (is.character(link)) 
      link <- match.arg(link, c("logit", "probit", "cloglog", "cauchit", "log", "loglog"))
    if (is.character(link.phi)) 
      link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))
    z$link <- list(mean = if (link != "loglog") make.link(link) else
                     structure(list(linkfun = function(mu) -log(-log(mu)), 
                                    linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 
                                      1 - .Machine$double.eps), .Machine$double.eps), 
                                    mu.eta = function(eta) {
                                      eta <- pmin(eta, 700)
                                      pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                                    }, valideta = function(eta) TRUE, name = "loglog"), 
                               class = "link-glm"),
                   precision = make.link(link.phi))
  }
  attr(z, "type") <- type
  z
}

FLXMRbeta_extra <- function(formula = .~., precision = ~ 1,
                      extra_components, offset = NULL,
                      link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                      link.phi = "log", control = betareg.control())
{
  if (!is(extra_components, "list")) extra_components <- list(extra_components)
  if (!all(sapply(extra_components, is, "extraComponent"))) stop("extra components need to be specified as a list of 'extraComponent' elements")
  new("FLXMRbeta_extra", FLXMRbeta(formula=formula, precision=precision, offset = offset,
                                link=link, link.phi=link.phi, control=control),
      extra_components = extra_components)
}

setMethod("FLXremoveComponent", signature(model = "FLXMRbeta_extra"),
          function(model, nok, ...) {
            if (any(seq_along(model@extra_components) %in% nok)) {
              drop <- which(seq_along(model@extra_components) %in% nok)
              model@extra_components <- model@extra_components[-drop]
              if (length(drop) == length(model@extra_components))
                model <- as(model, "FLXMRbeta")
            }
            model
          })

setMethod("FLXmstep", signature(model = "FLXMRbeta_extra"),
          function(model, weights, ...) {
  defineComponent <-
    list(uniform = expression({
      predict <- function(x, z, ...) 
        list(mean = x %*% coef, delta = delta)
      logLik <- function(x, y, z, ...) {
        pars <- predict(x, z, ...)
        lower <- with(pars, mean - delta)
        upper <- with(pars, mean + delta)
        - log(2 * pars$delta) + log(as.integer(lower < y & upper > y))
      }
      new("FLXcomponent", parameters = list(mean = coef, delta = delta), logLik = logLik, predict = predict, 
          df = df)
    }),
         betareg = FLXMRbeta()@defineComponent)
  extra_components <- lapply(model@extra_components,
                       function(x)
                             if (attr(x, "type") == "uniform")
                             with(x, with(list(coef = c(coef, rep(0, length.out = ncol(model@x) - length(coef))),
                                               df = 0, offset = NULL, delta = delta),
                                          eval(defineComponent[["uniform"]])))
                             else
                             with(x, with(list(coef = list(mean = c(coef$mean, rep(0, length.out = ncol(model@x) - length(coef$mean))),
                                                 precision = c(coef$precision, rep(0, length.out = ncol(model@z) - length(coef$precision)))),
                                               df = 0, offset = NULL, linkobjs = link),
                                          eval(defineComponent[["betareg"]]))))
  c(FLXmstep(as(model, "FLXMRbeta"), weights[, seq_len(ncol(weights) - length(model@extra_components)), drop=FALSE]),
    extra_components)
})

setMethod("FLXgetParameters", signature(object = "FLXMRbeta_extra"),
function(object, components) 
  callNextMethod(object, components[seq_len(length(components) - length(object@extra_components))]))

setMethod("FLXgetDesign", signature(object = "FLXMRbeta_extra"),
function(object, components) 
  FLXgetDesign(as(object, "FLXMRbeta"), components[seq_len(length(components) - length(object@extra_components))]))

setMethod("FLXreplaceParameters", signature(object="FLXMRbeta_extra"),
function(object, components, parms)
          c(FLXreplaceParameters(as(object, "FLXMRbeta"), components[seq_len(length(components) - length(object@extra_components))], parms),
            components[length(components) + seq(1 - length(object@extra_components), 0)]))

setMethod("FLXgradlogLikfun", signature(object="FLXMRbeta_extra"),
function(object, components, weights, ...)
          FLXgradlogLikfun(as(object, "FLXMRbeta"),
                           components[seq_len(length(components) - length(object@extra_components))],
                           weights[,seq_len(length(components) - length(object@extra_components)),drop=FALSE]))

setMethod("refit_optim", signature(object = "FLXMRbeta_extra"),
function(object, components, ...) {
  x <- refit_optim(as(object, "FLXMRbeta"), components[seq_len(length(components) - length(object@extra_components))], ...)
  names(x) <- paste("Comp", seq_along(x), sep = ".")
  x
})
