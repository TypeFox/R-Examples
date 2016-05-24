parallelSubscales <- function(dat, convertToNumeric = TRUE) {
  res <- list(input = list(dat = dat),
              intermediate = list(),
              output = list());
  
  if (!is.list(dat)) {
    stop("Argument 'dat' must be a dataframe or a list of dataframes!");
  }
  
  if (is.data.frame(dat)) {
    dat2 <- list();
    dat2[[deparse(substitute(dat))]] <- dat;
    dat <- dat2;
    rm(dat2);
  }

  res$intermediate$dat <- dat;
  
  ### Cycle through dataframes and average items over
  ### datapoints/observations/participants/records/rows
  ### and compute variances
  res$intermediate$itemMeans <- list();
  res$intermediate$itemVariances <- list();
  for (currentDf in 1:length(dat)) {
    res$intermediate$itemMeans[[currentDf]] <-
      colMeans(dat[[currentDf]], na.rm=TRUE);
    res$intermediate$itemVariances[[currentDf]] <-
      diag(var(dat[[currentDf]], na.rm=TRUE));
  }
  
  ### Make two dataframes out of these vectors
  res$intermediate$dat.itemMeans <- as.data.frame(res$intermediate$itemMeans);
  res$intermediate$dat.itemVariances <- as.data.frame(res$intermediate$itemVariances);
  
  ### Compute average over parallel tests
  res$intermediate$meanItemMeans <- rowMeans(res$intermediate$dat.itemMeans);
  res$intermediate$meanItemVariances <- rowMeans(res$intermediate$dat.itemVariances);

  ### Sort item means and variance
  res$intermediate$meanItemMeans.sorted <- sort(res$intermediate$meanItemMeans);
  res$intermediate$meanItemVariances.sorted <- sort(res$intermediate$meanItemVariances);
  
  ### Select items for test a1 and test a2
  res$output$byMeans <- list();
  res$output$byVariances <- list();
  
  res$output$byMeans$a1 <-
    names(res$intermediate$meanItemMeans.sorted)[is.odd(1:length(res$intermediate$meanItemMeans.sorted))];
  res$output$byMeans$a2 <-
    names(res$intermediate$meanItemMeans.sorted)[is.even(1:length(res$intermediate$meanItemMeans.sorted))];
  res$output$byVariances$a1 <-
    names(res$intermediate$meanItemVariances.sorted)[is.odd(1:length(res$intermediate$meanItemVariances.sorted))];
  res$output$byVariances$a2 <-
    names(res$intermediate$meanItemVariances.sorted)[is.even(1:length(res$intermediate$meanItemVariances.sorted))];
  
  ### Construct selection vectors
  res$output$byMeans$in_a1 <- names(res$intermediate$meanItemMeans) %in% res$output$byMeans$a1;
  res$output$byMeans$in_a2 <- names(res$intermediate$meanItemVariances) %in% res$output$byMeans$a2;
  res$output$byVariances$in_a1 <- names(res$intermediate$meanItemMeans) %in% res$output$byVariances$a1;
  res$output$byVariances$in_a2 <- names(res$intermediate$meanItemVariances) %in% res$output$byVariances$a2;
  
  ### Store means for each subscale
  res$output$byMeans$a1.mean <- mean(res$intermediate$meanItemMeans[res$output$byMeans$in_a1]);
  res$output$byMeans$a2.mean <- mean(res$intermediate$meanItemMeans[res$output$byMeans$in_a2]);
  res$output$byVariances$a1.mean <- mean(res$intermediate$meanItemMeans[res$output$byVariances$in_a1]);
  res$output$byVariances$a2.mean <- mean(res$intermediate$meanItemMeans[res$output$byVariances$in_a2]);

  ### Store variances for each subscale
  res$output$byMeans$a1.variance <- sum(res$intermediate$meanItemVariances[res$output$byMeans$in_a1]);
  res$output$byMeans$a2.variance <- sum(res$intermediate$meanItemVariances[res$output$byMeans$in_a2]);
  res$output$byVariances$a1.variance <- sum(res$intermediate$meanItemVariances[res$output$byVariances$in_a1]);
  res$output$byVariances$a2.variance <- sum(res$intermediate$meanItemVariances[res$output$byVariances$in_a2]);
  
  class(res) <- "parallelSubscales";
  return(res);

}

print.parallelSubscales <- function(x, nsmall=2, ...) {
    cat("--- When splitting by means:\n");
    cat(paste0("Mean subscale a1 = ",round(x$output$byMeans$a1.mean, digits=nsmall),
               ", variance = ", round(x$output$byMeans$a1.variance, digits=nsmall),
               " (", paste(x$output$byMeans$a1, collapse=", "), ")\n"));
    cat(paste0("Mean subscale a2 = ", round(x$output$byMeans$a2.mean, digits=nsmall),
               ", variance = ", round(x$output$byMeans$a2.variance, digits=nsmall),
               " (", paste(x$output$byMeans$a2, collapse=", "), ")\n"));
    cat("\n");
    cat("--- When splitting by variances:\n");
    cat(paste0("Mean subscale a1 = ", round(x$output$byVariances$a1.mean, digits=nsmall),
               ", variance = ", round(x$output$byVariances$a1.variance, digits=nsmall),
               " (", paste(x$output$byVariances$a1, collapse=", "), ")\n"));
    cat(paste0("Mean subscale a2 = ", round(x$output$byVariances$a2.mean, digits=nsmall),
               ", variance = ", round(x$output$byVariances$a2.variance, digits=nsmall),
               " (", paste(x$output$byVariances$a2, collapse=", "), ")\n"));
}


testRetestCES <- function(dat = NULL, moments = NULL,
                          testDat = NULL, retestDat = NULL,
                          parallelTests = 'means',
                          sortItems = FALSE, convertToNumeric = TRUE,
                          digits=4) {
  
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
  ### and retestDat with the data from the second administration. We now need
  ### to split this into the two parallel subscales. Ideally, the user provided
  ### a vector that distinguishes the items; alternatively, we will split the
  ### scale based on either the means or the variances.

  ### Regardless of whether it's necessary, we'll order this 'automated'
  ### generation of the parallel subscales.
  res$intermediate$parallelSubscales <-
    parallelSubscales(list(time1=res$intermediate$testDat,
                           time2=res$intermediate$retestDat));
  
  if (tolower(parallelTests) == "means") {
    res$intermediate$in_a1 <-
      res$intermediate$parallelSubscales$output$byMeans$in_a1;
    res$intermediate$in_a2 <-
      res$intermediate$parallelSubscales$output$byMeans$in_a2;
  }
  else if (tolower(parallelTests) == "variances") {
    res$intermediate$in_a1 <-
      res$intermediate$parallelSubscales$output$byVariances$in_a1;
    res$intermediate$in_a2 <-
      res$intermediate$parallelSubscales$output$byVariances$in_a2;
  }
  else if (length(parallelTests) > 1) {
    res$intermediate$in_a1 <- parallelTests == min(parallelTests);
    res$intermediate$in_a2 <- parallelTests != min(parallelTests);
  }
  else {
    stop("Argument 'parallelTests' must be either a vector identifying the two ",
         "parallel tests, or one of 'means' or 'variances'.");
  }
  
  ### Now we get the final four components (A at time1, A at time2,
  ### B at time1, and B at time2).
  res$intermediate$a1_time1.dat <- testDat[res$intermediate$in_a1];
  res$intermediate$a1_time2.dat <- retestDat[res$intermediate$in_a1];
  res$intermediate$a2_time1.dat <- testDat[res$intermediate$in_a2];
  res$intermediate$a2_time2.dat <- retestDat[res$intermediate$in_a2];
  
  ### Now we compute the actual scale scores from these separate items
  res$intermediate$a1_time1 <- rowSums(res$intermediate$a1_time1.dat);
  res$intermediate$a1_time2 <- rowSums(res$intermediate$a1_time2.dat);
  res$intermediate$a2_time1 <- rowSums(res$intermediate$a2_time1.dat);
  res$intermediate$a2_time2 <- rowSums(res$intermediate$a2_time2.dat);
  
  ### Compute the subscale means and variances to enable inspection of the
  ### assumption that the subscales are parallel.
  res$output$a1_time1_mean     <- mean(res$intermediate$a1_time1);
  res$output$a2_time1_mean     <- mean(res$intermediate$a2_time1);
  res$output$a1_time2_mean     <- mean(res$intermediate$a1_time2);
  res$output$a2_time2_mean     <- mean(res$intermediate$a2_time2);
  res$output$a1_time1_variance <- var(res$intermediate$a1_time1);
  res$output$a2_time1_variance <- var(res$intermediate$a2_time1);
  res$output$a1_time2_variance <- var(res$intermediate$a1_time2);
  res$output$a2_time2_variance <- var(res$intermediate$a2_time2);
  
  ### Store complete subscale variances; first the scale itself
  res$intermediate$A_time1 <- rowSums(testDat);
  res$intermediate$A_time2 <- rowSums(retestDat);
  ### The variance
  res$intermediate$A_time1.var <- var(res$intermediate$A_time1);
  res$intermediate$A_time2.var <- var(res$intermediate$A_time2);
  
  ### Store covariance and correlation matrices
  res$intermediate$cov <- cov(cbind(res$intermediate$a1_time1,
                                    res$intermediate$a1_time2,
                                    res$intermediate$a2_time1,
                                    res$intermediate$a2_time2));
  res$intermediate$cor <- cor(cbind(res$intermediate$a1_time1,
                                    res$intermediate$a1_time2,
                                    res$intermediate$a2_time1,
                                    res$intermediate$a2_time2));
  colnames(res$intermediate$cor) <- rownames(res$intermediate$cor) <-
    colnames(res$intermediate$cov) <- rownames(res$intermediate$cov) <-
    c('a1_time1', 'a1_time2', 'a2_time1', 'a2_time2');
  
  ### Check whether the two subscales have equal numbers of items
  if (ncol(res$intermediate$a1_time1.dat) == ncol(res$intermediate$a2_time1.dat)) {
    ### Equation 16 on page 212 of Schmidt, Le & Ilies (2003) defines the
    ### Coefficient of Equivalence and Stability (CES) as 2 times the sum of
    ### the covariance between a1_time1 & a1_time2 and the covariance between
    ### a1_time2 and a2_time1, divided by the product of the square roots of
    ### the variances of A_time1 and A_time2.
    
    res$output$testRetestCES <- (2 *
                                 (res$intermediate$cov['a1_time1', 'a2_time2'] +
                                  res$intermediate$cov['a1_time2', 'a2_time1'])) /
                                (sqrt(res$intermediate$A_time1.var) *
                                 sqrt(res$intermediate$A_time2.var));
    
    res$output$usedEquation <- "Equation 16 on page 212 of Schmidt, Le & Ilies (2003)"
  }
  else {
    ### Equation 17a on page 212 of Schmidt, Le & Ilies (2003) defines the
    ### Coefficient of Equivalence and Stability (CES) as the sum of the
    ### covariance between a1_time1 & a1_time2 and the covariance between
    ### a1_time2 and a2_time1, divided by the product of 2, p1, p2, the
    ### variance of sqrt(A_time1) and sqrt(A_time2), where p1 and p2 represent
    ### the proportions of items of the a1 and a2 subscales, respectively.
    
    res$intermediate$p1 <- p1 <- ncol(res$intermediate$a1_time1.dat) / ncol(testDat);
    res$intermediate$p2 <- p2 <-  ncol(res$intermediate$a2_time1.dat) / ncol(testDat);
    
    res$output$testRetestCES <- (res$intermediate$cov['a1_time1', 'a2_time2'] +
                                 res$intermediate$cov['a1_time2', 'a2_time1']) /
                                (2 * p1 * p2 * (sqrt(res$intermediate$A_time1.var) *
                                                sqrt(res$intermediate$A_time2.var)));
    
    res$output$usedEquation <- "Equation 17a on page 212 of Schmidt, Le & Ilies (2003)"
  }
  
  class(res) <- 'testRetestCES';
  return(res);
  
}

print.testRetestCES <- function(x, digits=x$input$digits, ...) {
  subscaleGeneration <- x$input$parallelTests;
  automatically <- "automatically ";
  if (length(subscaleGeneration) > 1) {
    subscaleGeneration <- "manually specified subscales";
    automatically <- "";
  }
  
  cat(paste0("Coefficient of Equivalence and Stability: ", round(x$output$testRetestCES, digits=digits), "\n\n",
             "To help assess whether the subscales (", automatically,
             "generated using ", subscaleGeneration,
             ") are parallel, here are the means and variances:\n",
             "Mean subscale a1, time 1: ", round(x$output$a1_time1_mean, digits=digits),
             " (variance = ", round(x$output$a1_time1_variance, digits=digits), ")\n",
             "Mean subscale a2, time 1: ", round(x$output$a2_time1_mean, digits=digits),
             " (variance = ", round(x$output$a2_time1_variance, digits=digits), ")\n",
             "Mean subscale a1, time 2: ", round(x$output$a1_time2_mean, digits=digits),
             " (variance = ", round(x$output$a1_time2_variance, digits=digits), ")\n",
             "Mean subscale a2, time 2: ", round(x$output$a2_time2_mean, digits=digits),
             " (variance = ", round(x$output$a2_time2_variance, digits=digits), ")\n"));
}
