### ROC S3 class from Thomas Lumley, RNews 2004.
### TestResult in every case is the statistic of interest, e.g. a t-statistic
### or a regression coefficient, and D is the 0-1 vector indicating
### whether it is a control (0) or a true finding (1).

print.ROC <- function(x,...){
  cat("ROC curve: ")
  print(x$call)
}

plot.ROC <- function(x, type = "b", null.line = TRUE,
                     xlab = "False Pos. Rate", ylab = "True Pos. Rate",
                     xlim = c(0, 1), ylim = c(0, 1), main = "ROC", ...)
{
  plot(x$mspec, x$sens, type = type, xlab = xlab, ylab = ylab,
       main = main, xlim = xlim, ylim = ylim, ...)
  if(null.line) abline(0, 1, lty = 3, col = "gray")

  invisible()
}

lines.ROC <- function(x,...){
  lines(x$mspec, x$sens, ...)
}

points.ROC <- function(x,...){
  points(x$mspec, x$sens, ...)
}

identify.ROC <- function(x, labels = NULL, ..., digits = 1)
{
  if (is.null(labels))
    labels <- round(x$test,digits)
  
  identify(x$mspec, x$sens, labels = labels,...)
}

ROC <- function(TestResult, ...) UseMethod("ROC")

ROC.default <- function(TestResult, D, take.abs = TRUE, ...){
  ## addition: D can also be a vector of indices
  if (length(D) < length(TestResult)) {
    D2 <- rep(0, length(TestResult))
    D2[D] <- 1
    D <- D2
  }

  if (take.abs) TestResult <- abs(TestResult)
  
  TT <- rev(sort(unique(TestResult)))
  DD <- table(-TestResult,D)
  
  sens <- cumsum(DD[,2])/sum(DD[,2])
  mspec <- cumsum(DD[,1])/sum(DD[,1])

  rval <- list(sens = sens, mspec = mspec,
               test = TT, call = sys.call())
  class(rval) <- "ROC"
  
  rval
}

AUC <- function(x, max.mspec = 1) {
  huhn <- aggregate(x$sens, list(x$mspec), max)
  mean(huhn[huhn[,1] <= max.mspec, 2])
}

roc.value <- function(found, true, totalN)
{
  TPR <- sum(found %in% true)
  FPR <- (length(found) - TPR) / (totalN - length(true))

  rval <- list(sens = TPR / length(true), mspec = FPR,
               test = NULL,
               call = sys.call())

  class(rval) <- "ROC"
  
  rval
}
