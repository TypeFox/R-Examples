lnl.negbin <- function(param, y, X, id, model, link, vlink, rn){

  K <- ncol(X)
  if (vlink == 'nb1') k <- 1 else k <- 2
  beta <- param[1L:K]
  if (model == "pooling") alpha <- param[K + 1L]
  if (model %in% c("between", "random")){
    a <- param[K + 1L]
    b <- param[K + 2L]
  }
  if (model == "within") a <- b <- 0
  bX <- as.numeric(crossprod(t(X), beta))
  l <- exp(bX)
  if (model == "pooling"){
    v <- l^(2-k)/alpha
    lnLi <- lgamma(y + v) - lgamma(y + 1) - lgamma(v) + v * log(v) -
      v * log(v + l) + y * log(l) - y * log(v + l)
  }
  else{
    Ti <- table(id)
    n <- length(unique(id))
    names.id <- as.character(unique(id))
    Li <- as.vector(tapply(l, id, sum))
    Yi <- as.vector(tapply(y, id, sum))
    names(Li) <- names(Yi) <- names.id
    lnA <- tapply(lgamma(l + y) - lgamma(l) - lgamma(y + 1), id, sum)
    lnB <- lgamma(Li) + lgamma(Yi + 1) - lgamma(Li + Yi)
    lnC <- lgamma(a + b) + lgamma(a + Li) + lgamma(b + Yi) -
      lgamma(a) - lgamma(b) - lgamma(a + b + Li + Yi)
    lnLi <- switch(model,
                   "within"  =   lnA + lnB,
                   "random"  =   lnA + lnC,
                   "between" = - lnB + lnC
                   )
  }
  lnL <- sum(lnLi)

  if (model == "pooling"){
    lb <- l
    vb <- (2 - k) * v
    vs <- - v / alpha
    Ll <- - (v + y) / (l + v) + y / l
    Lv <-  log(v) + 1 - log(l + v) - (v + y) / (l + v) + digamma(y + v) - digamma(v)
    gb <- (Ll * l + (2 - k) * v * Lv)
    gs <- - v / alpha * Lv
    gradi <-  cbind(gb * X, alpha = gs)
  }
  else{
    lnA.beta <- (digamma(l + y) - digamma(l)) * l
    lnB.beta <- (digamma(Li) - digamma(Li + Yi))[as.character(id)] * l
    lnC.beta <- (digamma(a + Li) - digamma(a + b + Li + Yi))[as.character(id)] * l
    lnC.a <- (digamma(a + b) + digamma(a + Li) - digamma(a) - digamma(a + b + Li + Yi)) / Ti
    lnC.b <- (digamma(a + b) + digamma(b + Yi) - digamma(b) - digamma(a + b + Li + Yi)) / Ti
    gradi <- switch(model,
                    "within"  =   lnA.beta + lnB.beta,
                    "random"  =   lnA.beta + lnC.beta,
                    "between" = - lnB.beta + lnC.beta
                    ) * X
    if (model != "within") gradi <- cbind(gradi, cbind(lnC.a, lnC.b)[as.character(id), ])
  }

  if (model == "pooling"){
    Lll <- (v + y) / (l + v)^2 - y / l^2
    Llv <- (y - l)/(l + v)^2
    Lvv <- 1 / v - 1 / (l + v) + (y - l) / (l + v)^2 + trigamma(y + v) - trigamma(v)
    lbb <- l
    vbb <- (2 - k)^2 * l
    vbs <- - (2 - k) * v / alpha
    vss <- 2 * v / alpha^2

    Hbb <- Ll * lbb + Lv * vbb +
      vb * (Lvv * vb + Llv * lb) +
        lb * (Llv * vb + Lll * lb)

    Hbv <- Lv * vbs +
      vb * (Lvv * vs ) +
        lb * (Llv * vs)

    Hvv <- Lv * vss +
      vs * (Lvv * vs)
    
    Hbb <- crossprod(X * Hbb, X)
    Hbv <- apply(Hbv * X, 2, sum)
    Hvv <- sum(Hvv)
    H <- rbind(cbind(Hbb, Hbv), c(Hbv, Hvv))
  }
  else{
    lXit <- l * X
    lXi <- apply(lXit, 2, tapply, id, sum)
    lnA.beta.beta <- (digamma(l + y) - digamma(l)) * l + 
      (trigamma(l + y) - trigamma(l)) * l^2
    lnA.beta.beta[is.na(lnA.beta.beta) | lnA.beta.beta<0] <- 0
    lnA.beta.beta <- crossprod(sqrt(lnA.beta.beta) * X)
    lnB.beta.beta.1 <- (digamma(Li) - digamma(Li + Yi))[as.character(id)] * l
    lnB.beta.beta.1 <-  - crossprod(sqrt( - lnB.beta.beta.1) * X)
    lnB.beta.beta.2 <- trigamma(Li) - trigamma(Li + Yi)
    lnB.beta.beta.2 <- crossprod(sqrt(lnB.beta.beta.2) * lXi)
    lnC.beta.beta.1 <- (digamma(a + Li) - digamma(a + b + Li + Yi))[as.character(id)] * l
    lnC.beta.beta.1 <-  - crossprod(sqrt( - lnC.beta.beta.1) * X)
    lnC.beta.beta.2 <- trigamma(a + Li) - trigamma(a + b + Li + Yi)
    lnC.beta.beta.2 <- crossprod(sqrt(lnC.beta.beta.2) * lXi)
    lnC.beta.a <- apply((trigamma(a + Li) - trigamma(a + b + Li + Yi)) * lXi,2,sum)
    lnC.beta.b <- apply( - trigamma(a + b + Li + Yi) * lXi,2,sum)
    lnC.a.a <- sum(trigamma(a + b) + trigamma(a + Li) - trigamma(a) - trigamma(a + b + Li + Yi))
    lnC.a.b <- sum(trigamma(a + b) - trigamma(a + b + Li + Yi))
    lnC.b.b <- sum(trigamma(a + b) + trigamma(b + Yi) - trigamma(b) - trigamma(a + b + Li + Yi))
    A.22 <- matrix(c(lnC.a.a,lnC.a.b,lnC.a.b,lnC.b.b),2,2)
    A.12 <- rbind(lnC.beta.a,lnC.beta.b)
    H <- switch(model,
                "within"=lnA.beta.beta + lnB.beta.beta.1 + lnB.beta.beta.2,
                "random"=lnA.beta.beta + lnC.beta.beta.1 + lnC.beta.beta.2,
                "between"= - lnB.beta.beta.1 - lnB.beta.beta.2 +
                lnC.beta.beta.1 + lnC.beta.beta.2
                )
    if (model!="within") H <- cbind(rbind(H, A.12), rbind(t(A.12), A.22))
  }
  attr(lnL, "gradient") <- gradi
  attr(lnL, "hessian") <-  H
  attr(lnL, "fitted.values") <- l
  lnL
}
