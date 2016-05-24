predict.ecoX <- function(object, newdraw = NULL, subset = NULL,
                         newdata = NULL, cond = FALSE, verbose = FALSE, ...){

  if (is.null(newdraw) && is.null(object$mu))
    stop("Posterior draws of mu and Sigma must be supplied")
  else if (!is.null(newdraw)){
    if (is.null(newdraw$mu) && is.null(newdraw$Sigma))
      stop("Posterior draws of both mu and Sigma must be supplied.")
    object <- newdraw
  }

  if (cond) { ## conditional prediction
    mu <- coef(object, subset = subset)
    n.draws <- nrow(mu)
    if (is.null(subset))
      subset <- 1:n.draws
    Sigma <- object$Sigma[subset,]
    if (is.null(newdata))
      X <- object$X
    else {
      mf <- match.call()
      if (is.matrix(eval.parent(mf$newdata)))
        data <- as.data.frame(data)
      tt <- terms(object)
      attr(tt, "intercept") <- 0
      X <- model.matrix(tt, newdata)
    }
    n <- nrow(X)
    res <- .C("preBaseX", as.double(X), as.double(mu), as.double(t(Sigma)),
              as.integer(length(c(X))), as.integer(nrow(mu)),
              as.integer(verbose),
              pdStore = double(n.draws*n*2), PACKAGE="eco")$pdStore
    res <- array(res, c(2, n, n.draws), dimnames=list(c("W1", "W2"),
                                          rownames(X), 1:n.draws))
    class(res) <- c("predict.ecoX", "array")
  }
  else {
    res <- predict.eco(object, newdraw = newdraw, subset = subset,
                       newdata = newdata, verbose = verbose, ...)
    colnames(res) <- c("W1", "W2", "X")
  }
  return(res)
}
