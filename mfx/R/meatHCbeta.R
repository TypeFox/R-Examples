meatHCbeta <-
function(x, type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", 
                                  "HC4", "HC4m", "HC5"), omega = NULL){
  X <- as.matrix(data.frame(x$x))
  if (any(alias <- is.na(coef(x)))) 
    X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)
  diaghat <- try(hatvalues(x), silent = TRUE)
  df <- n - NCOL(X)
  res <- rowMeans(estfun(x)/X, na.rm = TRUE)
  if (is.null(omega)) {
    type <- match.arg(type)
    if (type == "HC") 
      type <- "HC0"
    switch(type, const = {
      omega <- function(residuals, diaghat, df) rep(1, 
                                                    length(residuals)) * sum(residuals^2)/df
    }, HC0 = {
      omega <- function(residuals, diaghat, df) residuals^2
    }, HC1 = {
      omega <- function(residuals, diaghat, df) residuals^2 * 
        length(residuals)/df
    }, HC2 = {
      omega <- function(residuals, diaghat, df) residuals^2/(1 - 
                                                               diaghat)
    }, HC3 = {
      omega <- function(residuals, diaghat, df) residuals^2/(1 - 
                                                               diaghat)^2
    }, HC4 = {
      omega <- function(residuals, diaghat, df) {
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(4, n * diaghat/p)
        residuals^2/(1 - diaghat)^delta
      }
    }, HC4m = {
      omega <- function(residuals, diaghat, df) {
        gamma <- c(1, 1.5)
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], 
                                                      n * diaghat/p)
        residuals^2/(1 - diaghat)^delta
      }
    }, HC5 = {
      omega <- function(residuals, diaghat, df) {
        k <- 0.7
        n <- length(residuals)
        p <- as.integer(round(sum(diaghat), digits = 0))
        delta <- pmin(n * diaghat/p, pmax(4, n * k * 
                                            max(diaghat)/p))
        residuals^2/sqrt((1 - diaghat)^delta)
      }
    })
  }
  if (is.function(omega)) 
    omega <- omega(res, diaghat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  return(rval)
}
