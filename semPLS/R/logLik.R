### Armin (2012-09-12)

logLik.semplslv <- function (LV, object, REML = FALSE, ...) 
{
    res <- semPLS:::residuals.sempls(object)[, LV]
    pred <- semPLS:::predecessors(object$model)[[LV]]
    p <- length(pred)
    N <- object$N
    if (is.null(w <- object$obsweights)) {
        w <- rep.int(1, N)
    }
    else {
        excl <- w == 0
        if (any(excl)) {
            res <- res[!excl]
            N <- length(res)
            w <- w[!excl]
        }
    }
    N0 <- N
    if (REML) 
        N <- N - p
    val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + 
        log(sum(w * res^2))))
    if (REML){
        object$qr <- qr(object$factor_scores[, pred])
        val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
    }
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    val
}

logLik.sempls <- function(object, REML = FALSE, ..., LV){
  if(missing(LV)){
    endo <- semPLS:::endogenous(object$model)
    ret <- lapply(endo, logLik.semplslv, object, REML, ...)
    names(ret) <- endo
    class(ret) <- "logLiksempls"
    return(ret)
  }
  else{
    logLik.semplslv(LV, object, REML, ...)
  }
}

print.logLiksempls <- function(x, ...){
  out <- sapply(x, function(x) c(`log Lik.` = x, df = attr(x, "df")))
  out <- t(out)
  print(out, ...)
}

AIC.sempls <- function(..., k = 2, LV){
  objl <- list(...)
  stopifnot(all(sapply(objl, class) == "sempls"))
  ## names(objl) <- paste("AIC(m", seq_along(objl), ")", sep = "")
  names(objl) <- paste("AIC(", names(objl), ")", sep = "")
  if(missing(LV)){
    aic <- function(x){lapply(logLik(x), stats:::AIC.logLik, k = k)}
  }
  else{
    aic <- function(x){stats:::AIC.logLik(logLik(x, k = k, LV = LV))}
  }
  ll <- mapply(aic, objl)
  return(ll)
}


BIC.sempls <- function(..., k = 2, LV){
  objl <- list(...)
  stopifnot(all(sapply(objl, class) == "sempls"))
  ## names(objl) <- paste("BIC(m", seq_along(objl), ")", sep = "")
  names(objl) <- paste("BIC(", names(objl), ")", sep = "")
  if(missing(LV)){
    bic <- function(x){lapply(logLik(x), stats:::BIC.logLik)}
  }
  else{
    bic <- function(x){stats:::BIC.logLik(logLik(x, LV = LV))}
  }
  ll <- mapply(bic, objl)
  return(ll)
}

deviance.sempls <- function(object, LV){
  if(missing(LV))
    LV <- endo <- semPLS:::endogenous(object$model)
  colSums(semPLS:::residuals.sempls(object)[, LV, drop = FALSE]^2)
}

extractAIC.semplslv <- function (LV, fit, scale = 0, k = 2, ...){
    n <- fit$N
    pred <- semPLS:::predecessors(fit$model)[[LV]]
    p <- length(pred)
    edf <- p
    RSS <- deviance(fit, LV)
    dev <- if (scale > 0) 
        RSS/scale - n
    else n * log(RSS/n)
    c(df = edf, AIC = dev + k * edf)
}

extractAIC.sempls <- function (fit, scale = 0, k = 2, ..., LV){
  if(missing(LV)){
    endo <- semPLS:::endogenous(fit$model)
    ret <- t(sapply(endo, extractAIC.semplslv, fit, scale, k, ...))
    ## names(ret) <- endo
    colnames(ret) <- c("df", "AIC")
    return(ret)
  }
  else{
    extractAIC.semplslv(LV, fit, scale, k, ...)
  }
}
