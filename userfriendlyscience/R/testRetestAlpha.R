testRetestAlpha <- function(dat = NULL, moments = NULL,
                            testDat = NULL, retestDat = NULL,
                            sortItems = FALSE, convertToNumeric = TRUE) {
  
  res <- list(input = list(dat = dat,
                           moments = moments,
                           testDat = testDat,
                           retestDat = retestDat,
                           sortItems = sortItems,
                           convertToNumeric = convertToNumeric),
              intermediate = list(), output = list());
  
  ### If no dataframe was specified, load it from an SPSS file
  if (is.null(dat) && is.null(testDat) && is.null(retestDat)) {
    dat <- getData(errorMessage=paste0("No dataframe specified, and no valid datafile selected in ",
                                       "the dialog I then showed to allow selection of a dataset.",
                                       "Original error:\n\n[defaultErrorMessage]"),
                   use.value.labels=FALSE);
  }
  
  if (!is.null(dat)) {
    if (is.null(res$intermediate$moments)) {
      res$intermediate$moments <- rep(c(0,1), each=(ncol(dat))/2);
    }
    momentsBoolean <- (res$intermediate$moments == min(res$intermediate$moments));
    res$intermediate$testDat <- testDat <- dat[, momentsBoolean];
    res$intermediate$retestDat <- retestDat <- dat[, !momentsBoolean];
  }
  else if (xor(is.null(testDat), is.null(retestDat))) {
    stop("Provide both testDat and retestDat; or, if you have all scores in one ",
         "dataframe, provide it as 'dat' argument!");
  }
  
  if (sortItems) {
    res$intermediate$testDat <- testDat <- testDat[, order(names(testDat))];
    res$intermediate$retestDat <- retestDat <- retestDat[, order(names(retestDat))];
  }

  if (convertToNumeric) {
    res$intermediate$testDat <- testDat <- massConvertToNumeric(testDat);
    res$intermediate$retestDat <- retestDat <- massConvertToNumeric(retestDat);
  }
  
  ### So now we have a testDat and a retestDat, so we can get started.
  
  ### First get the sum of the covariances between
  ### the different-time/different-item covariances (equation 14, page 94
  ### of Green, 2003).
  res$intermediate$covar <- cov(testDat, retestDat);
  ### We have to remove the sum of the variances of course.
  res$intermediate$covar.sum <-
    sum(res$intermediate$covar) - sum(diag(res$intermediate$covar));
  
  ### Divide by (J(J-1)), where J is the number of items (equation 14)
  res$intermediate$J <- J <- ncol(testDat);
  res$intermediate$itemTrueScoreVariance <-
    res$intermediate$covar.sum / (J * (J - 1));
  
  ### Get item true score variance (equation 15, which is the numerator
  ### in the test-retest alpha)
  res$intermediate$testTrueScoreVariance <-
    J^2 * res$intermediate$itemTrueScoreVariance;
  
  ### Get common scale variance (equation 16, which is the denominator
  ### in the test-retest alpha)
  
  ### First compute the scales themselves
  res$intermediate$testScale <- rowSums(testDat);
  res$intermediate$retestScale <- rowSums(retestDat);
  
  ### The the product of the standard deviations
  res$intermediate$commonScaleVariance <-
    sd(res$intermediate$testScale) * sd(res$intermediate$retestScale);
  
  ### Then compute the test-retest alpha coefficient
  res$output$testRetestAlpha <- res$intermediate$testTrueScoreVariance /
    res$intermediate$commonScaleVariance;
  
  class(res) <- 'testRetestAlpha';
  return(res);
  
}

print.testRetestAlpha <- function(x, ...) {
  print(x$output$testRetestAlpha, ...);
}
