make.dmu.deta <- function(linkstr) switch(linkstr,
  "logit" = {
    logit_link <- make.link("logit")
    function(eta) logit_link$mu.eta(eta) * (1 - 2 * logit_link$linkinv(eta))
  },
  "probit" = function(eta) -eta * pmax(dnorm(eta), .Machine$double.eps),
  "cauchit" = function(eta) -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps),
  "cloglog" = function(eta) pmax((1 - exp(eta)) * exp(eta - exp(eta)), .Machine$double.eps),
  "loglog" = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
  "identity" = function(eta) rep.int(0, length(eta)),
  "log" = function(eta) pmax(exp(eta), .Machine$double.eps),
  "sqrt" = function(eta) rep.int(2, length(eta)),
  "1/mu^2" = function(eta) 3/(4 * eta^2.5),
  "inverse" = function(eta) 2/(eta^3)
)

