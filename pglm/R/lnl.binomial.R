lnl.binomial <- function(param, y, X, id, model, link, rn, start.sigma = FALSE){

  if (link == 'probit'){
    F <- pnorm
    f <- dnorm
    e <- function(x) - x * f(x)
  }
  if (link == 'logit'){
    F <- function(x) exp(x) / (1 + exp(x))
    f <- function(x) F(x) * (1 - F(x))
    e <- function(x) (exp(x) - exp(3*x))/(1+exp(x))^4
  }
  mills <- function(x) f(x) / F(x)
  millsp <- function(x) e(x)/F(x) - (f(x)/F(x))^2

  K <- ncol(X)
  n <- length(unique(id))
  N <- length(y)
  q <- 2 * y - 1

  beta <- param[1:K]
  bX <- as.numeric(crossprod(t(X), beta))
  if (start.sigma){
    mu <- -  tapply(q * mills(q * bX), id, sum) / tapply(millsp(q * bX), id, sum)
    return(sqrt(2) * sd(mu))
  }
  if (model == "random"){
    sigma <- param[K + 1L]
    Pitr <- lapply(rn$nodes, function(x) F( q * (bX + sigma * x)))
    Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
    Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE)) / sqrt(pi)
  }
  if (model == "pooling") Li <- F(q * bX)
  lnL <- sum(log(Li))

  if (model == "random"){
    gitr <- mapply(function(w, x, p)
                   q * w * cbind(rep(1, N), x) *
                   as.numeric(p[as.character(id)]) *
                   mills( q * (bX + sigma * x)),
                   rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
    g <- Reduce("+", gitr)
    gradi <- cbind(g[, 1] * X, g[, 2]) /
      as.numeric(Li[as.character(id)])/ sqrt(pi)
  }
  if (model == "pooling") gradi <- q * mills(q * bX) * X


  if (model == "pooling")
    H <- crossprod(millsp(q * bX) * X, X)
  if (model == "random"){
    Hr <- mapply(
                 function(w, v, p){
                   p <- p / Li
                   P <- p[as.character(id)]
                   p <- as.numeric(p)
                   P <- as.numeric(P)
                   z <- q * (bX + sigma * v)
                   gi <- q * mills(z) * cbind(X, v)
                   gi <- apply(gi, 2, tapply, id, sum)
                   H1 <- crossprod(p * gi, gi)
                   H2 <- crossprod(P * millsp(z) * cbind(X, v), cbind(X, v))
                   (H1 + H2) * w
                 },
                 rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
    H <- Reduce("+", Hr) / sqrt(pi) - crossprod(apply(gradi, 2, tapply, id, sum))
  }                  
  attr(lnL, 'gradient') <- gradi
  attr(lnL, 'hessian') <-   H
  lnL
}  


