#' Ordinal Probit Regression for Ordered Categorical Dependent Variables
#'
#' Vignette: \url{http://docs.zeligproject.org/en/latest/zeligchoice-oprobit.html}
#' @import methods
#' @export Zelig-ologit
#' @exportClass Zelig-ologit
#' 
#' @include model-obinchoice.R

zoprobit <- setRefClass("Zelig-oprobit",
                       contains = "Zelig-obinchoice")

zoprobit$methods(
  initialize = function() {
    callSuper()
    .self$name <- "oprobit"
    .self$packageauthors <- "William N. Venables, and Brian D. Ripley"
    .self$description <- "Ordinal Probit Regression for Ordered Categorical Dependent Variables"
    .self$method <- "probit"
    .self$linkinv <- function(eta, zeta) {
      tmp1 <- matrix(1, nrow = length(eta), ncol = length(zeta) + 1)
      tmp1[, 1:length(zeta)] <- pnorm(zeta - eta)
      return(tmp1)
    }
    .self$wrapper <- "oprobit"
    .self$vignette.url <- "http://docs.zeligproject.org/en/latest/zeligchoice-oprobit.html"
  }
)


zoprobit$methods(
  mcfun = function(x, b0=0, b1=1, ..., sim=TRUE){
    mu <- b0 + b1 * x
    n.sim = length(x)
    y.star <- rnorm(n = n.sim, mean = mu, sd = 1)  # latent continuous y
    t <- c(0,1,2)  # vector of cutpoints dividing latent space into ordered outcomes
    
    if(sim){
        y.obs <- rep(1, n.sim)
        for(i in 1:length(t)){
            y.obs <- y.obs + as.numeric(y.star > t[i]) # observed ordered outcome
        }
        return(as.factor(y.obs))
    }else{
        y.obs.hat <- rep(1, n.sim)
        for(i in 1:length(t)){
            y.obs.hat <- y.obs.hat + pnorm(q = t[i], mean = mu , sd = 1, lower.tail = FALSE) # expectation of observed ordered outcome
        }
        return(y.obs.hat)
    }
  }
)

