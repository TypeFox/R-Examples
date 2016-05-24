setClass("FLXMRlmc",
         representation(family = "character",
                        group = "factor",
                        censored = "formula",
                        C = "matrix"),
         contains = "FLXMR")

setClass("FLXMRlmcfix",
         contains = "FLXMRlmc")

setClass("FLXMRlmmc",
         representation(random = "formula",
                        z = "matrix",
                        which = "ANY"),
         contains = "FLXMRlmc")

setClass("FLXMRlmmcfix",
         contains = "FLXMRlmmc")

update.Residual <- function(fit, w, z, C, which, random, censored) {
  index <- lapply(C, function(x) x == 1)
  W <- rep(w, sapply(which, function(x) nrow(z[[x]])))
  ZGammaZ <- sapply(seq_along(which), function(i) sum(diag(crossprod(z[[which[i]]]) %*% random$Gamma[[i]])))
  WHICH <- which(sapply(C, sum) > 0)
  Residual <- if (length(WHICH) > 0)
    sum(sapply(WHICH, function(i) 
               w[i] * sum(diag(censored$Sigma[[i]]) - 2 * z[[which[i]]][index[[i]],,drop=FALSE] * censored$psi[[i]]))) else 0
  (sum(W * resid(fit)^2) + Residual + sum(w * ZGammaZ))/sum(W)
}

update.latent <- function(x, y, C, fit) {
  AnyMissing <- which(sapply(C, sum) > 0)
  index <- lapply(C, function(x) x == 1)
  Sig <- lapply(seq_along(x), function(i) fit$sigma2 * diag(nrow = nrow(x[[i]])))
  SIGMA <- rep(list(matrix(nrow = 0, ncol = 0)), length(x))
  if (length(AnyMissing) > 0) {
    SIGMA[AnyMissing] <- lapply(AnyMissing, function(i) {
      S <- Sig[[i]]
      SIG <- S[index[[i]], index[[i]]]
      if (sum(!index[[i]]) > 0) SIG <-
          SIG - S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*% S[!index[[i]],index[[i]]]
      SIG
    })
  }
  Sigma <- MU <- rep(list(vector("numeric", length = 0)), length(x))
  if (length(AnyMissing) > 0) {
    MU[AnyMissing] <- lapply(AnyMissing, function(i) {
      S <- Sig[[i]]
      Mu <- x[[i]][index[[i]],,drop=FALSE] %*% fit$coef
      if (sum(!index[[i]]) > 0) Mu <- Mu + S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*%
        (y[[i]][!index[[i]]] - x[[i]][!index[[i]],,drop=FALSE] %*% fit$coef)
      Mu
    })
  }
  moments <- lapply(seq_along(x), function(i) {
    if (sum(index[[i]]) > 0) moments_truncated(MU[[i]], SIGMA[[i]], y[[i]][C[[i]] == 1])
  })
  Sigma <- lapply(moments, "[[", "variance")
  censored <- list(mu = lapply(moments, "[[", "mean"),
                   Sigma = Sigma)
  list(censored = censored)
}

update.latent.random <- function(x, y, z, C, which, fit) {
  index <- lapply(C, function(x) x == 1)
  AnyMissing <- which(sapply(C, sum) > 0)
  Residual <- fit$sigma2$Residual
  Psi <- fit$sigma2$Random
  EVbeta <- lapply(seq_along(z), function(i) solve(1/Residual * crossprod(z[[i]]) + solve(Psi)))
  Sig <- lapply(seq_along(z), function(i) z[[i]] %*% Psi %*% t(z[[i]]) + Residual * diag(nrow = nrow(z[[i]])))
  SIGMA <- rep(list(matrix(nrow = 0, ncol = 0)), length(x))
  if (length(AnyMissing) > 0) {
    SIGMA[AnyMissing] <- lapply(AnyMissing, function(i) {
      S <- Sig[[which[i]]]
      SIG <- S[index[[i]], index[[i]]]
      if (sum(!index[[i]]) > 0) SIG <-
        SIG - S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*% S[!index[[i]],index[[i]]]
      SIG
    })
  }
  Sigma <- MU <- rep(list(vector("numeric", length = 0)), length(x))
  if (length(AnyMissing) > 0) {
    MU[AnyMissing] <- lapply(AnyMissing, function(i) {
      S <- Sig[[which[i]]]
      Mu <- x[[i]][index[[i]],,drop=FALSE] %*% fit$coef
      if (sum(!index[[i]]) > 0) {
        Mu <- Mu + S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*%
          (y[[i]][!index[[i]]] - x[[i]][!index[[i]],,drop=FALSE] %*% fit$coef)
      }
      Mu
    })
  }
  moments <- lapply(seq_along(x), function(i) {
    if (sum(index[[i]]) > 0) moments_truncated(MU[[i]], SIGMA[[i]], y[[i]][C[[i]] == 1])
  })
  Sigma <- lapply(moments, "[[", "variance")
  censored <- list(mu = lapply(moments, "[[", "mean"),
                   Sigma = Sigma,
                   psi = lapply(seq_along(x), function(i) {
                     if (sum(index[[i]]) > 0) return(diag(Sigma[[i]] %*% 
                                                          z[[which[i]]][index[[i]],,drop=FALSE] %*% EVbeta[[which[i]]])/Residual)
                     else return(vector("numeric", length = 0))
                   }))
  ybar <- lapply(seq_along(y), function(i) {
    Y <- y[[i]]
    Y[index[[i]]] <- censored$mu[[i]]
    Y
  })
  random <- list(beta = lapply(seq_along(x), function(i) 
                   EVbeta[[which[i]]] %*% t(z[[which[i]]]) %*% (ybar[[i]] - x[[i]] %*% fit$coef)/Residual),
                 Gamma = lapply(seq_along(x), function(i) {
                   if (sum(index[[i]]) > 0) {
                       return(EVbeta[[which[i]]] + (EVbeta[[which[i]]] %*% 
                                                    (t(z[[which[i]]][index[[i]],,drop=FALSE]) %*% censored$Sigma[[i]] %*%
                                                     z[[which[i]]][index[[i]],,drop=FALSE]) %*%
                                                    t(EVbeta[[which[i]]]))/Residual^2)
                     } else return(EVbeta[[which[i]]])
                 }))
  list(random = random,
       censored = censored)
}

moments_truncated <- function(mu, Sigma, T, ...) {
  Sigma <- as.matrix(Sigma)
  mu <- as.vector(mu)
  T <- as.vector(T)
  S <- 1/sqrt(diag(Sigma))
  T1 <- S * (T - mu)
  if (length(mu) == 1) {
    alpha <- pnorm(T1)
    dT1 <- dnorm(T1)
    Ex <- - dT1 / alpha
    Ex2 <- 1 - T1 * dT1 / alpha
  } else {
    R <- S * Sigma * rep(S, each = ncol(Sigma))
    diag(R) <- 1L
    alpha <- mvtnorm::pmvnorm(upper = T1, sigma = R, ...)
    rq <- lapply(seq_along(T1), function(q) (R - tcrossprod(R[,q])))
    R2 <- R^2
    Vq <- 1 - R2
    Sq <- sqrt(Vq)
    Rq <- lapply(seq_along(T1), function(q) rq[[q]]/(tcrossprod(Sq[,q])))
    Tq <- lapply(seq_along(T1), function(q) (T1 - R[,q] * T1[q])/Sq[,q])
    Phiq <- if (length(mu) == 1) 1 else
    sapply(seq_along(Rq), function(q) mvtnorm::pmvnorm(upper = Tq[[q]][-q], sigma = Rq[[q]][-q,-q], ...))
    phi_Phiq <- dnorm(T1) * Phiq
    Ex <- - (R %*% phi_Phiq)/alpha
    T2_entries <- lapply(seq_along(T1),
                         function(j) sapply(lapply(seq_along(T1)[seq_len(j)], function(i) R[,i] * T1 * phi_Phiq), function(z) sum(z * R[j,])))
    T2 <- diag(length(T1))
    T2[upper.tri(T2, diag = TRUE)] <- unlist(T2_entries)
    T2[lower.tri(T2)] <- t(T2)[lower.tri(T2)]
    phiqr <- lapply(seq_along(T1), function(q)
                    sapply(seq_along(T1), function(r) {
                      if (r == q) return(0) else return(mvtnorm::dmvnorm(T1[c(q, r)],
                            mean = rep(0, length.out = length(c(q,r))), sigma = R[c(q,r), c(q,r)]))}))
    if (length(mu) == 2) {
      Ex2 <- R - T2 / alpha +
        Reduce("+", lapply(seq_along(Tq), function(q)
                           tcrossprod(R[q,],
                                      rowSums(sapply(seq_along(Tq)[-q], function(r)
                                                     phiqr[[q]][r] * (R[,r] - R[q,r] * R[q,])))))) / alpha
    } else {
      betaq <- lapply(seq_along(T1), function(q) sweep(rq[[q]], 2, Vq[,q], "/"))
      Rqr <- lapply(seq_along(T1), function(q) 
                    lapply(seq_along(T1), function(r) if (r == q) return(0) else 
                           return((Rq[[q]][-c(q,r),-c(q,r)] - tcrossprod(Rq[[q]][-c(q,r),r]))/tcrossprod(sqrt(1 - Rq[[q]][-c(q,r),r]^2)))))
      Tqr <- lapply(seq_along(T1), function(q) {
        lapply(seq_along(T1), function(r) if (r == q) return(0) else 
               return((T1[-c(q,r)] - betaq[[r]][-c(r,q),q] * T1[q] - betaq[[q]][-c(r,q),r] * T1[r])/
                      (Sq[-c(r,q),q] * sqrt(1 - Rq[[q]][-c(r,q),r]^2))))})
      T3 <- Reduce("+", lapply(seq_along(Tq), function(q)
                               tcrossprod(R[q,],
                                          rowSums(sapply(seq_along(Tq)[-q], function(r)
                                                         phiqr[[q]][r] * (R[,r] - R[q,r] * R[q,]) *
                                                         mvtnorm::pmvnorm(upper = Tqr[[q]][[r]], sigma = Rqr[[q]][[r]], ...)))))) / alpha
      Ex2 <- R - T2 / alpha + 1/2 * (T3 + t(T3))

    }
  }
  moments <- list(mean = 1/S * Ex + mu,
                  variance = diag(1/S, nrow = length(T)) %*% (Ex2 - tcrossprod(Ex)) %*% diag(1/S, nrow = length(T)))
  if (!all(is.finite(unlist(moments))) || any(moments$mean > T) || any(eigen(moments$variance)$values < 0)) {
      moments <- list(mean = T - abs(diag(Sigma)),
                      variance = Sigma)
  }
  moments
}

FLXMRlmmc <- function(formula = . ~ ., random, censored, varFix, eps = 10^-6, ...)
{
  family <- "gaussian"
  censored <- if (length(censored) == 3) censored else formula(paste(".", paste(deparse(censored), collapse = "")))
  if (missing(random)) {
    if (missing(varFix)) varFix <- FALSE
    else if ((length(varFix) > 1) || (is.na(as.logical(varFix)))) stop("varFix has to be a logical vector of length one")
    object <- new("FLXMRlmc", formula = formula, censored = censored,
                  weighted = TRUE, family = family, name = "FLXMRlmc:gaussian")
    if (varFix) object <- new("FLXMRlmcfix", object)
    lmc.wfit <- function(x, y, w, C, censored) {
      W <- rep(w, sapply(x, nrow))
      X <- do.call("rbind", x)
      AnyMissing <- which(sapply(C, sum) > 0)
      ybar <- lapply(seq_along(y), function(i) {
        Y <- y[[i]]
        Y[C[[i]] == 1] <- censored$mu[[i]]
        Y
      })
      Y <- do.call("rbind", ybar)
      fit <- lm.wfit(X, Y, W, ...)
      fit$sigma2 <- if (length(AnyMissing) > 0) (sum(W * resid(fit)^2) +
                                                 sum(sapply(AnyMissing, function(i) 
                                                            w[i] * sum(diag(censored$Sigma[[i]])))))/sum(W) else sum(W * resid(fit)^2)/sum(W)
      fit$df <- ncol(X)
      fit
    }
    object@defineComponent <- expression({
      predict <- function(x, ...)
        lapply(x, function(X) X %*% coef)
      
      logLik <- function(x, y, C, group, censored, ...) {
        AnyMissing <- which(sapply(C, sum) > 0)
        index <- lapply(C, function(x) x == 1)
        V <- lapply(x, function(X) diag(nrow(X)) * sigma2)
        mu <- predict(x, ...)
        SIGMA <- rep(list(matrix(nrow = 0, ncol = 0)), length(x))
        if (length(AnyMissing) > 0) {
          SIGMA[AnyMissing] <- lapply(AnyMissing, function(i) {
            S <- V[[i]]
            SIG <- S[index[[i]], index[[i]]]
            if (sum(!index[[i]]) > 0) SIG <-
              SIG - S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*% S[!index[[i]],index[[i]]]
            SIG
          })
        }
        MU <- rep(list(vector("numeric", length = 0)), length(x))
        if (length(AnyMissing) > 0) {
          MU[AnyMissing] <- lapply(AnyMissing, function(i) {
            S <- V[[i]]
            Mu <- mu[[i]][index[[i]]]
            if (sum(!index[[i]]) > 0) Mu <- Mu + S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*%
              (y[[i]][!index[[i]]] - mu[[i]][!index[[i]]])
            Mu
          })
        }
        llh <- sapply(seq_along(x), function(i) {
          LLH <- 0
          if (sum(index[[i]]) > 0) LLH <- log(mvtnorm::pmvnorm(upper = y[[i]][index[[i]]], mean = as.vector(MU[[i]]),
                                                               sigma = SIGMA[[i]]))
          if (sum(!index[[i]]) > 0) LLH <- LLH + mvtnorm::dmvnorm(t(y[[i]][!index[[i]]]), mean = mu[[i]][!index[[i]]],
                                                                  sigma = V[[i]][!index[[i]], !index[[i]], drop = FALSE], log=TRUE)
          LLH/nrow(V[[i]])
        })
        as.vector(ungroupPriors(matrix(llh), group, !duplicated(group)))
      }
      new("FLXcomponent", parameters = list(coef = coef, sigma2 = sigma2, 
            censored = censored), logLik = logLik, predict = predict, 
            df = df)
    })
    object@fit <- if (varFix) {
      function(x, y, w, C, fit) {
        any_removed <- any(w <= eps)
        if (any_removed) {
          ok <- apply(w, 2, function(x) x > eps)
          W <- lapply(seq_len(ncol(ok)), function(i) w[ok[,i],i])
          X <- lapply(seq_len(ncol(ok)), function(i) x[ok[,i],,drop = FALSE])
          y <- lapply(seq_len(ncol(ok)), function(i) y[ok[,i]])
          C <- lapply(seq_len(ncol(ok)), function(i) C[ok[,i]])
        } else {
          X <- rep(list(x), ncol(w))
          y <- rep(list(y), ncol(w))
          C <- rep(list(C), ncol(w))
          W <- lapply(seq_len(ncol(w)), function(i) w[,i])
        }
        if ("coef" %in% names(fit[[1]]))
          fit <- lapply(seq_len(ncol(w)), function(k) update.latent(X[[k]], y[[k]], C[[k]], fit[[k]]))
        else {
          fit <- lapply(seq_len(ncol(w)), function(k)
                        list(censored = list(mu = lapply(seq_along(y[[k]]), function(i) y[[k]][[i]][C[[k]][[i]] == 1]),
                                 Sigma = lapply(C[[k]], function(x) diag(1, nrow = sum(x)) * var(unlist(y[[k]]))))))
        }
        fit <- lapply(seq_len(ncol(w)), function(k) c(lmc.wfit(X[[k]], y[[k]], W[[k]], C[[k]], fit[[k]]$censored),
                                                      censored = list(fit[[k]]$censored)))
        sigma2 <- sum(sapply(fit, function(x) x$sigma2) * colMeans(w))
        for (k in seq_len(ncol(w))) fit[[k]]$sigma2 <- sigma2
        lapply(fit, function(Z) with(list(coef = coef(Z),
                                          df = Z$df + 1/ncol(w),
                                          sigma2 =  Z$sigma2,
                                          censored = Z$censored),
                                     eval(object@defineComponent)))
      }
    } else {
      function(x, y, w, C, fit){
        any_removed <- any(w <= eps)
        if (any_removed) {
          ok <- w > eps
          w <- w[ok]
          x <- x[ok,,drop = FALSE]
          y <- y[ok]
          C <- C[ok]
        }
        if ("coef" %in% names(fit)) {
          fit <- update.latent(x, y, C, fit)
        } else {
          fit$censored <- list(mu = lapply(seq_along(y), function(i) y[[i]][C[[i]] == 1]),
                               Sigma = lapply(C, function(x) diag(1, nrow = sum(x)) * var(unlist(y))))
        }
        fit <- c(lmc.wfit(x, y, w, C, fit$censored),
                 censored = list(fit$censored))
        with(list(coef = coef(fit),
                  df = fit$df + 1,
                  sigma2 =  fit$sigma2,
                  censored = fit$censored),
             eval(object@defineComponent))
      }
    }
  } else {
    if (missing(varFix)) varFix <- c(Random = FALSE, Residual = FALSE)
    else if (length(varFix) != 2 || is.null(names(varFix)) || any(is.na(pmatch(names(varFix), c("Random", "Residual"))))) 
      stop("varFix has to be a named vector of length two")
    else names(varFix) <- c("Random", "Residual")[pmatch(names(varFix), c("Random", "Residual"))]
    random <- if (length(random) == 3) random else formula(paste(".", paste(deparse(random), collapse = "")))
    object <- new("FLXMRlmmc", formula = formula, random = random, censored = censored,
                  weighted = TRUE, family = family, name = "FLXMRlmmc:gaussian")
    if (any(varFix)) object <- new("FLXMRlmmcfix", object)
    add <- function(x) Reduce("+", x)    
    lmmc.wfit <- function(x, y, w, z, C, which, random, censored) {
      effect <- lapply(seq_along(which), function(i) z[[which[i]]] %*% random$beta[[i]])
      Effect <- do.call("rbind", effect)
      W <- rep(w, sapply(x, nrow))
      X <- do.call("rbind", x)
      ybar <- lapply(seq_along(y), function(i) {
        Y <- y[[i]]
        Y[C[[i]] == 1] <- censored$mu[[i]]
        Y
      })
      Y <- do.call("rbind", ybar)
      fit <- lm.wfit(X, Y - Effect, W, ...)
      wGamma <- add(lapply(seq_along(which), function(i) w[i] * random$Gamma[[i]]))
      bb <- add(lapply(seq_along(which), function(i) tcrossprod(random$beta[[i]]) * w[i]))
      fit$sigma2 <- list(Random = (wGamma + bb)/sum(w))
      fit$df <- ncol(X)
      fit
    }
    
    object@defineComponent <- expression({
      predict <- function(x, ...)
        lapply(x, function(X) X %*% coef)
      
      logLik <- function(x, y, z, C, which, group, censored, ...) {
        AnyMissing <- which(sapply(C, sum) > 0)
        index <- lapply(C, function(x) x == 1)
        V <- lapply(z, function(Z) tcrossprod(tcrossprod(Z, sigma2$Random), Z) + diag(nrow(Z)) * sigma2$Residual)
        mu <- predict(x, ...)
        SIGMA <- rep(list(matrix(nrow = 0, ncol = 0)), length(x))
        if (length(AnyMissing) > 0) {
          SIGMA[AnyMissing] <- lapply(AnyMissing, function(i) {
            S <- V[[which[i]]]
            SIG <- S[index[[i]], index[[i]]]
            if (sum(!index[[i]]) > 0) SIG <-
              SIG - S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*% S[!index[[i]],index[[i]]]
            SIG
          })
        }
        MU <- rep(list(vector("numeric", length = 0)), length(x))
        if (length(AnyMissing) > 0) {
          MU[AnyMissing] <- lapply(AnyMissing, function(i) {
            S <- V[[which[i]]]
            Mu <- mu[[i]][index[[i]]]
            if (sum(!index[[i]]) > 0) Mu <- Mu + S[index[[i]],!index[[i]]] %*% solve(S[!index[[i]],!index[[i]]]) %*%
              (y[[i]][!index[[i]]] - mu[[i]][!index[[i]]])
            Mu
          })
        }
        llh <- sapply(seq_along(x), function(i) {
          LLH <- 0
          if (sum(index[[i]]) > 0) LLH <- log(mvtnorm::pmvnorm(upper = y[[i]][index[[i]]], mean = as.vector(MU[[i]]),
                                                               sigma = SIGMA[[i]]))
          if (sum(!index[[i]]) > 0) LLH <- LLH + mvtnorm::dmvnorm(t(y[[i]][!index[[i]]]), mean = mu[[i]][!index[[i]]],
                                                                sigma = V[[which[i]]][!index[[i]], !index[[i]], drop = FALSE], log=TRUE)
          LLH/nrow(V[[which[i]]])
      })
        as.vector(ungroupPriors(matrix(llh), group, !duplicated(group)))
      }
      new("FLXcomponent", parameters = list(coef = coef, sigma2 = sigma2, 
            censored = censored, random = random), logLik = logLik, predict = predict, 
            df = df)
    })
    object@fit <- if (any(varFix)) {
      function(x, y, w, z, C, which, fit) {
        any_removed <- any(w <= eps)
        if (any_removed) {
          ok <- apply(w, 2, function(x) x > eps)
          W <- lapply(seq_len(ncol(ok)), function(i) w[ok[,i],i])
          X <- lapply(seq_len(ncol(ok)), function(i) x[ok[,i],,drop = FALSE])
          y <- lapply(seq_len(ncol(ok)), function(i) y[ok[,i]])
          C <- lapply(seq_len(ncol(ok)), function(i) C[ok[,i]])
          which <- lapply(seq_len(ncol(ok)), function(i) which[ok[,i]])
        } else {
          X <- rep(list(x), ncol(w))
          y <- rep(list(y), ncol(w))
          C <- rep(list(C), ncol(w))
          which <- rep(list(which), ncol(w))
          W <- lapply(seq_len(ncol(w)), function(i) w[,i])
        }
        if ("coef" %in% names(fit[[1]])) 
          fit <- lapply(seq_len(ncol(w)), function(k) update.latent.random(X[[k]], y[[k]], z, C[[k]], which[[k]],
                                                                           fit[[k]]))
        else {
            fit <- lapply(seq_len(ncol(w)), function(k)
                          list(random = list(beta = lapply(W[[k]], function(i) rep(0, ncol(z[[i]]))),
                                              Gamma = lapply(W[[k]], function(i) diag(ncol(z[[i]])))),
                               censored = list(mu = lapply(seq_along(y[[k]]), function(i) y[[k]][[i]][C[[k]][[i]] == 1]),
                                                Sigma = lapply(C[[k]], function(x) diag(1, nrow = sum(x)) * var(unlist(y[[k]]))),
                                                psi = lapply(C[[k]], function(x) rep(0, sum(x))))))
        }
        fit <- lapply(seq_len(ncol(w)), function(k) c(lmmc.wfit(X[[k]], y[[k]], W[[k]], z, C[[k]],
                                                                which[[k]], fit[[k]]$random, fit[[k]]$censored),
                                                      random = list(fit[[k]]$random),
                                                      censored = list(fit[[k]]$censored)))
        if (varFix["Random"]) {
          prior_w <- apply(w, 2, weighted.mean, w = sapply(x, length))
          Psi <- add(lapply(seq_len(ncol(w)), function(k) fit[[k]]$sigma2$Random * prior_w[k]))
          for (k in seq_len(ncol(w))) fit[[k]]$sigma2$Random <- Psi
        }
        for (k in seq_len(ncol(w))) 
          fit[[k]]$sigma2$Residual <- update.Residual(fit[[k]], W[[k]], z, C[[k]], which[[k]],
                                                      fit[[k]]$random, fit[[k]]$censored)
        if (varFix["Residual"]) {
          prior <- colMeans(w)
          Residual <- sum(sapply(fit[[k]]$sigma2$Residual, function(x) x) * prior)
          for (k in seq_len(ncol(w))) fit[[k]]$sigma2$Residual <- Residual
        } 
        n <- nrow(fit[[1]]$sigma2$Random)
        lapply(fit, function(Z) with(list(coef = coef(Z),
                                          df = Z$df + n*(n+1)/(2*ifelse(varFix["Random"], ncol(w), 1)) +
                                          ifelse(varFix["Residual"], 1/ncol(w), 1),
                                          sigma2 =  Z$sigma2,
                                          random = Z$random,
                                          censored = Z$censored),
                                     eval(object@defineComponent)))
      }
    } else {
      function(x, y, w, z, C, which, fit){
        any_removed <- any(w <= eps)
        if (any_removed) {
          ok <- w > eps
          w <- w[ok]
          x <- x[ok,,drop = FALSE]
          y <- y[ok]
          C <- C[ok]
          which <- which[ok]
        }
        if ("coef" %in% names(fit)) fit <- update.latent.random(x, y, z, C, which, fit)
        else {
            fit <- list(random = list(beta = lapply(which, function(i) rep(0, ncol(z[[i]]))),
                            Gamma = lapply(which, function(i) diag(ncol(z[[i]])))),
                        censored = list(mu = lapply(seq_along(y), function(i) y[[i]][C[[i]] == 1]),
                            Sigma = lapply(C, function(x) diag(1, nrow = sum(x)) * var(unlist(y))),
                            psi = lapply(C, function(x) rep(0, sum(x)))))
        }
        fit <- c(lmmc.wfit(x, y, w, z, C, which, fit$random, fit$censored),
                 random = list(fit$random), censored = list(fit$censored))
        fit$sigma2$Residual <- update.Residual(fit, w, z, C, which, fit$random, fit$censored)
        n <- nrow(fit$sigma2$Random)
        with(list(coef = coef(fit),
                  df = fit$df + n*(n+1)/2 + 1,
                  sigma2 =  fit$sigma2,
                  random = fit$random,
                  censored = fit$censored),
             eval(object@defineComponent))
      }
    }
  }
  object
}

setMethod("FLXmstep", signature(model = "FLXMRlmc"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  return(sapply(1:ncol(weights), function(k) model@fit(model@x, model@y, weights[,k], model@C,  
                                                       components[[k]]@parameters)))
})

setMethod("FLXmstep", signature(model = "FLXMRlmcfix"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  return(model@fit(model@x, model@y, weights, model@C, 
                   lapply(components, function(x) x@parameters)))
})

setMethod("FLXmstep", signature(model = "FLXMRlmmc"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  return(sapply(1:ncol(weights), function(k) model@fit(model@x, model@y, weights[,k], model@z, model@C, model@which, 
                                                       components[[k]]@parameters)))
})

setMethod("FLXmstep", signature(model = "FLXMRlmmcfix"),
          function(model, weights, components)
{
  weights <- weights[!duplicated(model@group),,drop=FALSE]
  return(model@fit(model@x, model@y, weights, model@z, model@C, model@which,
                   lapply(components, function(x) x@parameters)))
})

setMethod("FLXgetModelmatrix", signature(model="FLXMRlmc"),
          function(model, data, formula, lhs=TRUE, ...)
{
  formula_nogrouping <- RemoveGrouping(formula)
  if (formula_nogrouping == formula) stop("please specify a grouping variable")
  model <- callNextMethod(model, data, formula, lhs)
  model@fullformula <- update(model@fullformula,
                              paste(".~. |", .FLXgetGroupingVar(formula)))
  mt2 <- terms(model@censored, data=data)
  mf2 <- model.frame(delete.response(mt2), data=data, na.action = NULL)
  model@C <- model.matrix(attr(mf2, "terms"), data)
  model@group <- grouping <- .FLXgetGrouping(formula, data)$group
  model@x <- matrix(lapply(unique(grouping), function(g) model@x[grouping == g, , drop = FALSE]), ncol = 1)
  if (lhs) model@y <- matrix(lapply(unique(grouping), function(g) model@y[grouping == g, , drop = FALSE]), ncol = 1)
  model@C <- matrix(lapply(unique(grouping), function(g) model@C[grouping == g, , drop = FALSE]), ncol = 1)
  model
})

setMethod("FLXgetModelmatrix", signature(model="FLXMRlmmc"),
          function(model, data, formula, lhs=TRUE, ...)
{
  model <- callNextMethod(model, data, formula, lhs)
  mt1 <- terms(model@random, data=data)
  mf1 <- model.frame(delete.response(mt1), data=data, na.action = NULL)
  model@z <- model.matrix(attr(mf1, "terms"), data)
  grouping <- .FLXgetGrouping(formula, data)$group
  z <- matrix(lapply(unique(grouping), function(g) model@z[grouping == g, , drop = FALSE]), ncol = 1)
  model@z <- unique(z)
  model@which <- match(z, model@z)
  model
})

setMethod("FLXgetObs", "FLXMRlmc", function(model) sum(sapply(model@x, nrow)))

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRlmc"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@C, model@group, x@parameters$censored))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRlmmc"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@z, model@C, model@which, model@group, x@parameters$censored))
})

setMethod("predict", signature(object="FLXMRlmc"), function(object, newdata, components, ...)
{
  object <- FLXgetModelmatrix(object, newdata, formula = object@fullformula, lhs = FALSE)
  lapply(components, function(comp) unlist(comp@predict(object@x, ...)))
})

