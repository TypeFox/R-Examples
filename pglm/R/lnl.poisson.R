lnl.poisson <- function(param, y, X, id, model, link, rn, other = NULL){

  if (is.null(other)) other <- "inv"
  
  Ti <- table(id)
  n <- length(unique(id))
  names.id <- as.character(unique(id))
  K <- ncol(X)
  beta <- param[1L:K]
  if (model != "pooling"){
    if (other != "sd") d <- param[K+1L]
    else d <- 1 / param[K+1L]
  }
  bX <- as.numeric(crossprod(t(X), beta))
  lit <- exp(bX)
  if (model == "pooling") lnLi <- - lit + y * bX - lgamma(y + 1)
  else{
    Li <- as.vector(tapply(lit, id, sum))
    Yi <- as.vector(tapply(y, id, sum))
    names(Li) <- names(Yi) <- names.id
    lnA <- tapply(y * log(lit) - lgamma(y + 1), id, sum)
    lnB <- lgamma(Yi + 1) - log(Li) * Yi
    lnC <- d * log(d) - (d + Yi) * log(Li + d) + lgamma(d + Yi) - lgamma(d)
    lnLi <- switch(model,
                   "within"  =   lnA + lnB,
                   "random"  =   lnA + lnC,
                   "between" = - lnB + lnC
                   )
  }
  lnL <- sum(lnLi)
  
  if (model == "pooling"){
    gradi <- (y - lit) * X
  }
  else{
    lnA.beta <- y
    lnB.beta <- - (Yi / Li)[as.character(id)] * lit
    lnC.beta <- - ((d + Yi)/(d + Li))[as.character(id)] * lit
    lnC.d <- (1 + log(d) - log(d + Li) - (d + Yi)/(d + Li) +
              digamma(d + Yi) - digamma(d)) / Ti
    if (other == "sd") lnC.d <- - d^2 * lnC.d
    gradi <- switch(model,
                    "within"  =   lnA.beta + lnB.beta,
                    "random"  =   lnA.beta + lnC.beta,
                    "between" = - lnB.beta + lnC.beta
                    )
    gradi <- gradi * X
    if (model != "within") gradi <- cbind(gradi, lnC.d[as.character(id)])
  }
  attr(lnL, 'gradient') <- gradi
  
  if (model == "pooling"){
    H <- crossprod( - lit * X, X)
  }
  else{
    lXit <- lit * X
    lXi <- apply(lXit, 2, tapply, id, sum)
    lnA.beta.beta <- lnA.beta.d <- 0
    lnB.beta.beta.1 <- -(Yi / Li)[as.character(id)]*lit
    lnB.beta.beta.2 <- (Yi / Li^2)
    lnB.beta.beta.1 <- -crossprod(sqrt(- lnB.beta.beta.1) * X)
    lnB.beta.beta.2 <- crossprod(sqrt(lnB.beta.beta.2) * lXi)
    lnB.beta.beta <- lnB.beta.beta.1 + lnB.beta.beta.2
    lnB.beta.d <- 0
    lnC.beta.beta.1 <- -((d+Yi)/(d+Li))[as.character(id)]*lit
    lnC.beta.beta.2 <- ((d+Yi)/(d+Li)^2)
    lnC.beta.beta.1 <- -crossprod(sqrt(-lnC.beta.beta.1)*X)
    lnC.beta.beta.2 <- crossprod(sqrt(lnC.beta.beta.2)*lXi)
    lnC.beta.beta <- lnC.beta.beta.1+lnC.beta.beta.2
    lnC.beta.d <- -((Li-Yi)/(d+Li)^2)[as.character(id)]*lit
    lnC.beta.d <- apply(lnC.beta.d*X,2,sum)  
    lnC.d.d <- (1/d-1/(d+Li)-(Li-Yi)/(d+Li)^2+trigamma(d+Yi)-trigamma(d))
    lnC.d.d <- sum(lnC.d.d)
    A.22 <- lnC.d.d
    A.12 <- lnC.beta.d
    if (other == "sd"){
      A.22 <- d^4 * A.22 - d * sum(lnC.d)
      A.12 <- - d^2 * A.12
    }
    H <- switch(model,
                "within"  = lnA.beta.beta + lnB.beta.beta,
                "random"  = lnA.beta.beta + lnC.beta.beta,
                "between" =-lnB.beta.beta + lnC.beta.beta
                )
    if (model != "within"){
      H <- cbind(
                 rbind(H, A.12),
                 c(t(A.12), A.22)
                 )
    }
  }
  attr(lnL, "hessian") <- H
  attr(lnL, "fitted.values") <- lit
  lnL
}
