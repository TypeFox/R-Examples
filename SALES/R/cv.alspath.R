cv.alspath <- function(outlist, lambda, x, y, foldid, 
  pred.loss, tau) {
  typenames <- "Asymmetric squared error loss"
  if (!match(pred.loss, c("loss"), FALSE)) {
    warning("Only 'loss' available for ALS regression; 'loss' used")
    pred.loss <- "loss"
  }
  y <- as.double(y)
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- ercls(y-predmat, tau)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
      1))
  list(cvm = cvm, cvsd = cvsd, name = typenames)
} 
