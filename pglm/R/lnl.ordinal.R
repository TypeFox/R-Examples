lnl.ordinal <- function(param, y, X, id, model, link, rn, start.sigma = FALSE){

  if (start.sigma) gradient <- hessian <- TRUE

  # F, f, e la fonction de distribution de proba et ses derivees
  # premiere et seconde
  if (link == "probit"){
    F <- pnorm
    f <- dnorm
    e <- function(x) - x * f(x)
  }
  if (link == "logit"){
    F <- function(x) exp(x)/(1+exp(x))
    f <- function(x) exp(x)/(1+exp(x))^2
    e <- function(x) (exp(x) - exp(3*x))/(1+exp(x))^4
  }
  K <- ncol(X)
  n <- length(unique(id))
  J <- length(unique(y))
  beta <- param[1L:K]
  mu <- c(-100, 0, param[(K+1L):(K+J-2L)], 100)
  bX <- as.numeric(crossprod(t(X), beta))
  if (model == "pooling"){
    z1 <- mu[y+1] - bX                 ; z2 <- mu[y] - bX
    F1 <- F(z1)                        ; F2 <- F(z2)
    Li <- F1 - F2
  }
  if (model == "random"){
    sigma <- param[K+J-1L]
    Pitr <- lapply(rn$nodes, function(x)
                   F(mu[y+1] - bX - sigma * x) - F(mu[y] - bX - sigma * x))
    Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
    Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE))/sqrt(pi)
  }
  lnL <- sum(log(Li))
  W <- matrix(0, nrow(X), J)
  # The first two cols are not estimated coefficients (-infty and 0)
  # so remove them
  W1 <- (col(W) == (y + 1))[, -c(1, 2)]  ; W2 <- (col(W) == y)[, -c(1, 2)]
  if (model == "pooling"){
    f1 <- f(mu[y+1] - bX)                ; f2 <- f(mu[y] - bX)
    gmu <- (W1 * f1 - W2 * f2)
    gb <- - (f1 - f2) * X
    gradi <- cbind(gb, gmu) / Li
  }
  if (model == "random"){
    gitr <- mapply(
                   function(w, x, p){
                     f1 <- f(mu[y+1] - bX - sigma * x) ; f2 <- f(mu[y] - bX - sigma * x)
                     F1 <- F(mu[y+1] - bX - sigma * x) ; F2 <- F(mu[y] - bX - sigma * x)
                     w * as.numeric(p[as.character(id)]) *
                       (f1 * cbind(- 1, W1, - x) -
                        f2 * cbind(- 1, W2, - x)) / (F1 - F2)
                   },
                   rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
    g <- Reduce("+", gitr)
    gradi <- cbind(X * g[, 1], g[, -c(1)])/as.numeric(Li[as.character(id)])/sqrt(pi)
  }

  if (model == "pooling"){
    e1 <- e(mu[y+1] - bX)              ; e2 <- e(mu[y] - bX)                     
    M1 <- cbind(- X, W1)
    M2 <- cbind(- X, W2)
    H <- crossprod(M1 * e1 / Li, M1) - crossprod(M2 * e2 / Li, M2) - crossprod(gradi)
  }
  if (model == "random"){
    H <- mapply(
                function(w, x, p){
                  p <- p / Li
                  P <- p[as.character(id)]
                  p <- as.numeric(p)
                  P <- as.numeric(P)
                  z1 <- mu[y+1] - bX - sigma * x
                  z2 <- mu[y] - bX - sigma * x
                  f1 <- f(z1) ; f2 <- f(z2)
                  e1 <- e(z1) ; e2 <- e(z2)
                  F1 <- F(z1) ; F2 <- F(z2)
                  M1 <- cbind(- X, W1, - x) ; M2 <- cbind(- X, W2, - x)
                  git <- (f1 * M1 - f2 * M2) / (F1 - F2)
                  gi <- apply(git, 2, tapply, id, sum)
                  H1 <- crossprod(p * gi, gi)
                  H2 <- crossprod(P * (e1 * M1) / (F1 - F2), M1) -
                    crossprod(P * (e2 * M2) / (F1 - F2), M2) -
                      crossprod( git * sqrt(P))
                  (H1 + H2) * w},
                rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
    H <- Reduce("+", H) / sqrt(pi) - crossprod(apply(gradi, 2, tapply, id, sum))
  }
  if (start.sigma){
    gradit <- (- f1 + f2 ) / Li
    hessit <- (e1 - e2) / Li - gradit^2
    mu <- tapply(gradit, id, sum) / tapply(hessit, id, sum)
    return(sqrt(2) * sd(mu))
  }
  attr(lnL, "gradient") <- gradi
  attr(lnL, "hessian") <- H
  lnL
}    
