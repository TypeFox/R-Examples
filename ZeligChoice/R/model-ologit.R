#' Ordinal Logistic Regression for Ordered Categorical Dependent Variables
#'
#' Vignette: \url{http://docs.zeligproject.org/en/latest/zeligchoice-ologit.html}
#' @import methods
#' @export Zelig-ologit
#' @exportClass Zelig-ologit
#' 
#' @include model-obinchoice.R

zologit <- setRefClass("Zelig-ologit",
                       contains = "Zelig-obinchoice")

zologit$methods(
  initialize = function() {
    callSuper()
    .self$name <- "ologit"
    .self$packageauthors <- "William N. Venables, and Brian D. Ripley"
    .self$description <- "Ordinal Logit Regression for Ordered Categorical Dependent Variables"
    .self$method <- "logistic"
    .self$linkinv <- function(eta, zeta) {
      tmp1 <- matrix(1, nrow = length(eta), ncol = length(zeta) + 1)
      tmp1[, 1:length(zeta)] <- exp(zeta - eta) / (1 + exp(zeta - eta))
      return(tmp1)
    }
    .self$wrapper <- "ologit"
    .self$vignette.url <- "http://docs.zeligproject.org/en/latest/zeligchoice-ologit.html"
  }
)


zologit$methods(
  mcfun = function(x, b0=0, b1=1, ..., sim=TRUE){
    mu <- b0 + b1 * x
    n.sim = length(x)
    y.star <- rlogis(n = n.sim, location = mu, scale = 1)  # latent continuous y 
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
        y.obs.hat <- y.obs.hat + plogis(q = t[i], location = mu , scale = 1, lower.tail = FALSE) # expectation of observed ordered outcome
      }
      return(y.obs.hat)
    }
  }
)

