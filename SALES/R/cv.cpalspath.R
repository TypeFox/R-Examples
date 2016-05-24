cv.cpalspath <- function(outlist, lambda, x, y, foldid, 
  pred.loss, w, tau) {
  typenames <- "Coupled asymmetric squared error loss"
  if (!match(pred.loss, c("loss"), FALSE)) {
    warning("Only 'loss' available for coupled ALS regression; 'loss' used")
    pred.loss <- "loss"
  }
  y <- as.double(y)
  nfolds <- max(foldid)
  predmat1 <- matrix(NA, length(y), length(lambda))
  predmat2 <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- (foldid == i)
    fitobj <- outlist[[i]]
    nlami <- length(fitobj$lambda)
    preds <- predict(fitobj, x[whichfold, , drop = FALSE], type = "response")
    predmat1[whichfold, seq(nlami)] <- preds
    predmat2[whichfold, seq(nlami)] <- predict(fitobj, x[whichfold, , drop = FALSE], 
      type = "scale") + preds
    nlams[i] <- nlami
  }
  cvraw <- w * ercls(y-predmat1, 0.5) + ercls(y-predmat2, tau)
  N <- length(y) - apply(is.na(predmat2), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
      1))
  list(cvm = cvm, cvsd = cvsd, name = typenames)
} 
