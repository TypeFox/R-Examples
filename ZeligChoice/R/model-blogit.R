#' Bivariate Logistic Regression for Two Dichotomous Dependent Variables
#'
#' Vignette: \url{http://docs.zeligproject.org/en/latest/zeligchoice-blogit.html}
#' @import methods
#' @export Zelig-blogit
#' @exportClass Zelig-blogit
#' 
#' @include model-bbinchoice.R

zblogit <- setRefClass("Zelig-blogit",
                       contains = "Zelig-bbinchoice")

zblogit$methods(
  initialize = function() {
    callSuper()
    .self$name <- "blogit"
    .self$description <- "Bivariate Logit Regression for Dichotomous Dependent Variables"
    .self$family <- quote(binom2.or(zero = 3))
    .self$linkinv <- binom2.or()@linkinv
    .self$wrapper <- "blogit"
    .self$vignette.url <- "http://docs.zeligproject.org/en/latest/zeligchoice-blogit.html"
  }
)

zblogit$methods(
  mcfun = function(x, b0=0, b1=1, b2=1, b3=0.5, ..., sim=TRUE){
    n.sim = length(x)
    pi1 <- 1/(1 + exp(b0 + b1 * x))
    pi2 <- 1/(1 + exp(b2 + b3 * x))

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