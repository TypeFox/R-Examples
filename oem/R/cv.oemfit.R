cv.oemfit <- function(formula, data = list(), lambda = NULL,
                      type.measure = c('mse', 'mae'), ...,
                      nfolds = 10, foldid,
                      penalty = c("lasso", "scad", "elastic.net",
                        "ngarrote", "mcp")) {
  this.call <- match.call()
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  N <- nrow(x)
  y <- drop(y) # only vector form is allowed
  if (missing(type.measure)) type.measure = "mse"
  else type.measure <- match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2)
    stop ("Need more than on value of lambda for cv.oemfit")
  oemfit.object <- oemfit(formula, data = data, lambda = lambda, penalty = penalty)
  lambda <- oemfit.object$lambda
  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = N))
  else nfolds <- max(foldid)
  nz <- sapply(predict(oemfit.object, type = 'nonzero'), length)
  if (nfolds < 3) stop("nfolds must be greater than 3; nfolds = 10 is recommended")
  outlist <- as.list(seq(nfolds))
  ######################3
  # fit the n-fold model and store
  for (i in seq(nfolds)) {
    index <- foldid == i
    y_sub <- y[!index]
    xx <- x[!index,, drop = FALSE]
    outlist[[i]] <- oemfit(y_sub ~ xx - 1, lambda = lambda, penalty = penalty)
  }
  # Use type.measure to evaluate
  typenames <- c(mse = "Mean-Squared Error", mae = "Mean Absolute Error")
  if (!match(type.measure, c("mse", "mae"), FALSE)){
    warning ("Only 'mse' and 'mae' available; 'mse' used")
    type.measure <- "mse"
  }
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    index <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[index,, drop = FALSE])
    nlami <- length(outlist[[i]]$lambda)
    predmat[index, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  N.cv <- length(y) - apply(is.na(predmat), 2, sum)
  cvraw <- switch(type.measure,
                  'mse' = (y - predmat)^2,
                  'mae' = abs(y - predmat)
                  )

  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)
               /(N.cv - 1))
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvup = cvm + cvsd, cvlo = cvm - cvsd,
              nzero = nz, fit = oemfit.object,
              name = typenames[type.measure])
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  obj$penalty <- penalty
  class(obj) <- "cv.oemfit"
  obj
}
