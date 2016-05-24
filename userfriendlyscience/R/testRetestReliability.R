testRetestReliability <- function(dat = NULL, moments = NULL,
                                  testDat = NULL, retestDat = NULL,
                                  parallelTests = 'means',
                                  sortItems = FALSE, convertToNumeric = TRUE,
                                  digits=2) {

  ### Make object to store results
  res <- list(input = list(dat = dat,
                           moments = moments,
                           testDat = testDat,
                           retestDat = retestDat,
                           parallelTests = parallelTests,
                           sortItems = sortItems,
                           convertToNumeric = convertToNumeric,
                           digits = digits),
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
      if (is.odd(ncol(dat))) {
        stop("If argument 'dat' is provided, this dataframe is split into two. ",
             "Therefore, it must have an even number of columns.");
      }
      res$intermediate$moments <- rep(c(0,1), each=(ncol(dat))/2);
    }
    res$intermediate$momentsBoolean <- (res$intermediate$moments == min(res$intermediate$moments));
    res$intermediate$testDat <- testDat <- dat[, res$intermediate$momentsBoolean];
    res$intermediate$retestDat <- retestDat <- dat[, !res$intermediate$momentsBoolean];
  }
  else if (xor(is.null(testDat), is.null(retestDat))) {
    stop("Provide both testDat and retestDat; or, if you have all scores in one ",
         "dataframe, provide it as 'dat' argument!");
  }

  ### Store number of observations
  res$input$n.observations <- nrow(testDat);
  
  if (sortItems) {
    res$intermediate$testDat <- testDat <- testDat[, order(names(testDat))];
    res$intermediate$retestDat <- retestDat <- retestDat[, order(names(retestDat))];
  }
  
  if (ncol(testDat) != ncol(retestDat)) {
    stop("The dataframe for each measurement moment must have the same number of ",
         "items. The current version of testRetestCES only supports compuring the ",
         "test-retest CES for a scale that is split into parallel halves post-hoc; ",
         "see Schmidt, Le & Ilies (2003), pages 210 and 212.");
  }
  
  if (ncol(testDat) < 2) {
    stop("The scale at each measurement moment must contain at least two ",
         "items to split into subscales. The scale you specified has only ",
         ncol(testDat), " items.");
  }
  
  if (convertToNumeric) {
    res$intermediate$testDat <- testDat <- massConvertToNumeric(testDat);
    res$intermediate$retestDat <- retestDat <- massConvertToNumeric(retestDat);
  }
  
  ### So, now we have testDat with the data from the first administration,
  ### and retestDat with the data from the second administration. We can
  ### now call testRetestAlpha and testRetestCES.
  
  res$intermediate$testRetestAlpha <-
    testRetestAlpha(testDat = testDat, retestDat = retestDat);
  res$intermediate$testRetestCES <-
    testRetestCES(testDat = testDat, retestDat = retestDat,
                  parallelTests = parallelTests);
  
  ### Extract test-retest alpha and CES and store it in the output object
  res$output$testRetestAlpha <- res$intermediate$testRetestAlpha$output$testRetestAlpha;
  res$output$testRetestCES <- res$intermediate$testRetestCES$output$testRetestCES;
  
  class(res) <- "testRetestReliability";
  ### Return result
  return(res);
}

print.testRetestReliability <- function (x, digits=x$input$digits, ...) {
  cat(paste0("                         Items at time 1: ", paste(names(x$intermediate$testDat), collapse=", "),
           "\n                         Items at time 2: ", paste(names(x$intermediate$retestDat), collapse=", "),
           "\n                            Observations: ", x$input$n.observations,
           "\n           Test-retest Alpha Coefficient: ", round(x$output$testRetestAlpha, digits=digits),
           "\n"));
  print(x$intermediate$testRetestCES, digits=digits);
  invisible();
}
