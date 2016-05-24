"bickde" <-
function(data, factor=10) {
  require(KernSmooth)
  hpre <- dpik(apply(data, 1, mean))
  cat("Preliminary choice (based on midpoints & DPI bandwidth):\n")
  print(hpre)
  cat("Likelihood Cross-Validation Bandwidth Choice for ICE Kernel DE\n \n")
  cat("Caution: This optimization can be very SLOW .... \n ")
  cat("Type CTRL-C if you need to escape!\n")
  h1 <- hpre/factor
  h2 <- hpre*factor
  cvout <- optimize(likelihoodcv, c(h1, h2),tol=h1, data=data)
  c(h = cvout$minimum)
}

