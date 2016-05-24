#' Bivariate Probit Regression for Two Dichotomous Dependent Variables
#'
#' Vignette: \url{http://docs.zeligproject.org/en/latest/zeligchoice-bprobit.html}
#' @import methods
#' @export Zelig-bprobit
#' @exportClass Zelig-bprobit
#' 
#' @include model-bbinchoice.R

zbprobit <- setRefClass("Zelig-bprobit",
                        contains = "Zelig-bbinchoice")

zbprobit$methods(
  initialize = function() {
    callSuper()
    .self$name <- "bprobit"
    .self$description <- "Bivariate Probit Regression for Dichotomous Dependent Variables"
    .self$family <- quote(binom2.rho(zero = 3))
    .self$linkinv <- binom2.rho()@linkinv
    .self$wrapper <- "bprobit"
    .self$vignette.url <- "http://docs.zeligproject.org/en/latest/zeligchoice-bprobit.html"
  }
)

zbprobit$methods(
  mcfun = function(x, b0=0, b1=1, b2=1, b3=0.5, ..., sim=TRUE){
    n.sim = length(x)
    pi1 <- pnorm(b0 + b1 * x)
    pi2 <- pnorm(b2 + b3 * x)

    if(sim){
      y1 <- rbinom(n=n.sim, size=1, prob=pi1)
      y2 <- rbinom(n=n.sim, size=1, prob=pi2)
      return(as.data.frame(y1, y2, x))
    }else{
      y1.hat <- pi1
      y2.hat <- pi2
      return(as.data.frame(y1.hat, y2.hat, x))
    }
  }
)