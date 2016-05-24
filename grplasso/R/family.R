setClass("grpl.model", representation = representation(
    invlink          = "function",
    link             = "function",                     
    nloglik          = "function",                   
    ngradient        = "function",
    nhessian         = "function",
    check            = "function",
    name             = "character",
    comment          = "character"
))

setMethod("show", "grpl.model", function(object) {
    cat("Model:", object@name, "\n")
    cat("Comment:", object@comment, "\n\n")
})

grpl.model <- function(invlink, link, nloglik, ngradient, nhessian,
                       check, name = "user-specified",
                       comment = "user-specified"){
  ## Purpose: Generates models to be used for the Group Lasso algorithm.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## invlink: a function with arguments "eta" implementing the inverse link
  ##          function.
  ## link:    a function with arguments "mu" implementing the link
  ##          function.
  ## nloglik: a function with arguments "y", "mu" and "weights"
  ##          implementing the negative log-likelihood function.
  ## ngradient: a function with arguments "x", "y", "mu" and "weights" 
  ##            implementing the negative gradient of the log-likelihood
  ##            function.
  ## nhessian: a function with arguments "x", "mu" and "weights"
  ##           implementing the negative} hessian of the log-likelihood
  ##           function.
  ## name: a character name
  ## comment: a character comment
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:12

    RET <- new("grpl.model",
               invlink   = invlink,
               link      = link,
               nloglik   = nloglik,
               ngradient = ngradient,
               nhessian  = nhessian,
               check     = check,
               name      = name,
               comment   = comment)
    RET
}

## Logistic Regression 
LogReg <- function(){
  grpl.model(invlink   = function(eta) 1 / (1 + exp(-eta)),
             link      = function(mu) log(mu / (1 - mu)),
             nloglik   = function(y, eta, weights, ...)
               -sum(weights * (y * eta - log(1 + exp(eta)))),
             ngradient = function(x, y, mu, weights, ...)
               -crossprod(x, weights * (y - mu)),
             nhessian  = function(x, mu, weights, ...)
               crossprod(x, weights * mu * (1 - mu) * x),
             check     = function(y) all(y %in% c(0, 1)),
             name      = "Logistic Regression Model",
             comment   = "Binary response y has to be encoded as 0 and 1")
}

## Linear Regression
LinReg <- function(){
  grpl.model(invlink  = function(eta) eta,
             link  = function(mu) mu,
             nloglik   = function(y, eta, weights, ...)
               sum(weights * (y - eta)^2),
             ngradient = function(x, y, mu, weights, ...)
               -2 * crossprod(x, weights * (y - mu)),
             nhessian  = function(x, mu, weights, ...)
               2 * crossprod(x, weights * x),
             check     = function(y) TRUE,
             name      = "Linear Regression Model",
             comment   = "Use update.hess=\"lambda\" in grpl.control because the Hessian is constant")
}

## Poisson Regression
PoissReg <- function(){
  grpl.model(invlink    = function(eta) exp(eta),
             link       = function(mu) log(mu),
             nloglik    = function(y, eta, weights, ...)
               sum(weights * (exp(eta) - y * eta)),
             ngradient  = function(x, y, mu, weights, ...)
               -crossprod(x, weights * (y - mu)),
             nhessian   = function(x, mu, weights, ...)
               crossprod(x, weights * mu * x),
             check      = function(y) all(y >= 0) & all(y == ceiling(y)),
             name       = "Poisson Regression Model",
             comment    = "")
}

