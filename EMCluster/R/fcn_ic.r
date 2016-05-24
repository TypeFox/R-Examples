### This file contains tool functions for information criteria.
### Written: Wei-Chen Chen on 2008/11/02.


em.ic <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL,
    llhdval = NULL){
  if(is.null(emobj)){
    emobj <- list(n = nrow(x), pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  emobj$llhdval <- logL(x, emobj = emobj)
  emobj$M <- length(emobj$pi) - 1 + length(emobj$Mu) + length(emobj$LTSigma)

  ret <- list()
  ret$AIC <- em.aic(x, emobj = emobj)
  ret$BIC <- em.bic(x, emobj = emobj)
  ret$ICL <- em.icl(x, emobj = emobj)
  ret$ICL.BIC <- em.icl.bic(x, emobj = emobj)
  ret$CLC <- em.clc(x, emobj = emobj)

  ret
}

em.aic <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  if(is.null(emobj$M)){
    emobj$M <- length(emobj$pi) - 1 + length(emobj$Mu) + length(emobj$LTSigma)
  }
  if(is.null(emobj$llhdval)){
    emobj$llhdval <- logL(x, emobj = emobj)
  }

  -2 * emobj$llhdval + 2 * emobj$M
}

em.bic <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  if(is.null(emobj$n)){
    emobj$n <- nrow(x)
  }
  if(is.null(emobj$M)){
    emobj$M <- length(emobj$pi) - 1 + length(emobj$Mu) + length(emobj$LTSigma)
  }
  if(is.null(emobj$llhdval)){
    emobj$llhdval <- logL(x, emobj = emobj)
  }

  -2 * emobj$llhdval + log(emobj$n) * emobj$M
}

em.icl <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  if(is.null(emobj$n)){
    emobj$n <- nrow(x)
  }
  if(is.null(emobj$M)){
    emobj$M <- length(emobj$pi) - 1 + length(emobj$Mu) + length(emobj$LTSigma)
  }
  if(is.null(emobj$llhdval)){
    emobj$llhdval <- logL(x, emobj = emobj)
  }

  Z.unnorm <- e.step(x, emobj = emobj, norm = FALSE)$Gamma
  logL.map <- do.call("c", lapply(1:emobj$n, function(i){ max(Z.unnorm[i, ]) }))

  -2 * sum(logL.map) + log(emobj$n) * emobj$M
}

em.icl.bic <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  if(is.null(emobj$n)){
    emobj$n <- nrow(x)
  }
  if(is.null(emobj$M)){
    emobj$M <- length(emobj$pi) - 1 + length(emobj$Mu) + length(emobj$LTSigma)
  }
  if(is.null(emobj$llhdval)){
    emobj$llhdval <- logL(x, emobj = emobj)
  }

  Z <- e.step(x, emobj = emobj)$Gamma
  Z.unnorm <- e.step(x, emobj = emobj, norm = FALSE)$Gamma
  log.Z <- Z.unnorm - log(dmixmvn(x, emobj = emobj))
  logL.EN <- Z * log.Z

  -2 * (emobj$llhdval + sum(logL.EN)) + log(emobj$n) * emobj$M
}

em.clc <- function(x, emobj = NULL, pi = NULL, Mu = NULL, LTSigma = NULL){
  if(is.null(emobj)){
    emobj <- list(pi = pi, Mu = Mu, LTSigma = LTSigma)
  }
  if(is.null(emobj$llhdval)){
    emobj$llhdval <- logL(x, emobj = emobj)
  }

  Z <- e.step(x, emobj = emobj)$Gamma
  Z.unnorm <- e.step(x, emobj = emobj, norm = FALSE)$Gamma
  log.Z <- Z.unnorm - log(dmixmvn(x, emobj = emobj))
  logL.EN <- Z * log.Z

  -2 * (emobj$llhdval + sum(logL.EN))
}

