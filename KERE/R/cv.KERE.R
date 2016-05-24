cv.KERE <- function(x, y, kern, lambda = NULL,
                    nfolds = 5, foldid, omega = 0.5, ...) {
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  y <- as.double(y)
  x <- as.matrix(x)
  N <- NROW(x)
  if (missing(foldid))
    foldid <-
    sample(rep(seq(nfolds), length = N))
  else
    nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=5 recommended")
  ###Now fit the nfold models and store them
  nlams <- double(nfolds)
  predmat <- matrix(NA, length(y), length(lambda))
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    fit <- KERE(
      x = x[!which, , drop = FALSE],
      y = y_sub, kern = kern, lambda = lambda, omega = omega, ...
    )
    preds <-
      predict.KERE(fit, kern, x[!which, ,drop = FALSE], x[which, , drop = FALSE])
    nlami <- length(fit$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- ercls(y - predmat, omega)
  cvob <- cvcompute(cvraw, foldid, nlams)
  cvraw <- cvob$cvraw
  N <- cvob$N
  cvm <- colMeans(cvraw, na.rm = TRUE)
  cvsd <-
    sqrt(colMeans(scale(cvraw, cvm, FALSE) ^ 2, na.rm = TRUE) / (N - 1))
  cvm <- na.omit(cvm)
  cvsd <- na.omit(cvsd)
  cvm.min <- min(cvm)
  idmin <- cvm <= cvm.min
  lambda.min <- max(lambda[idmin])
  obj <-
    list(
      lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm +
        cvsd, cvlo = cvm - cvsd, name = "Expectile Loss", cvm.min = cvm.min,
      lambda.min = lambda.min
    )
  class(obj) <- "cv.KERE"
  obj
}
