
samplingDistribution <- function(popValues = c(0, 1), popFrequencies = c(50, 50),
                                 sampleSize = NULL, sampleFromPop = FALSE, ...) {
  
  if (is.null(sampleSize)) {
    sampleSize <- sum(popFrequencies);
  }

  if (sampleFromPop) {
    sampleVector <- sample(popValues, size=sampleSize,
                           replace=TRUE, prob=popFrequencies);
  }
  else {
    sampleVector <- rep(popValues, times=popFrequencies);    
  }
  
  return(normalityAssessment(sampleVector = sampleVector, ...));
  
}

dataShape <- function(sampleVector, na.rm=TRUE, type=2, digits=2,
                      conf.level=.95) {
  
  ### Note: this is adapted from the 'skewness' and 'kurtosis' functions
  ### in 'e1071' and the way they're computed in 'describe' in 'psych', and
  ### for the variances, also based on Stan Browns page at
  ### http://www.tc3.edu/instruct/sbrown/stat/shape.htm, as well as of course
  ### on Joanes & Gill (1998)
  
  res <- list(input = list(sampleVector = sampleVector, na.rm=na.rm,
                           type=type, digits=digits),
              intermediate = list(), output = list());
  
  if (any(isNA <- is.na(sampleVector))) {
    if (na.rm) {
      res$intermediate$sampleVector <- sampleVector <- sampleVector[!isNA];
    }
    else {
      return(NA);
    }
  }
  
  if (!(type %in% (1:3))) {
    stop("Invalid 'type' argument. This must be 1, 2 or 3; see ?skewness for more info.");
  }
  
  ### Number of observations
  res$intermediate$n <- n <- length(sampleVector);
  
  if (n < 3) {
    stop("There are only ", n, " datapoints/observations - this is too few to ",
         "compute a sensible estimate of the 'shape' of this data.");
  }
  
  ### First moment (mean)
  res$intermediate$m1 <-       sum(sampleVector^1) / n;

  ### Center the data
  res$intermediate$centeredVector <- sampleVector <- sampleVector - res$intermediate$m1;

  ### Second moment (variance)
  res$intermediate$m2 <- m2 <- sum(sampleVector^2) / n;
  
  ### Third moment
  res$intermediate$m3 <- m3 <- sum(sampleVector^3) / n;
  
  ### Fourth moment
  res$intermediate$m4 <- m4 <- sum(sampleVector^4) / n;

  ### Sample skewness
  res$intermediate$g1 <- m3 / (m2 ^ (3/2));

  ### Sample kurtosis
  res$intermediate$g2 <- (m4 / (m2 ^ 2)) - 3;
  
  ###  Population skewness
  res$intermediate$G1 <- res$intermediate$g1 *
                         sqrt( n * (n-1) ) /
                         ( n - 2 );

  ### Population kurtosis
  res$intermediate$G2 <- (n-1) *
                         ( (n+1) * res$intermediate$g2 + 6 ) /
                         ( (n-2) * (n-3) );
  
  ### Type 3 correction for skewness
  res$intermediate$b1 <- res$intermediate$g1 *
                         (((n-1) / n) ^ (3/2));
  
  ### Type 3 correction for kurtosis
  res$intermediate$b2 <- (res$intermediate$g2 + 3) *
                         (((n-1) / n) ^ 2) - 3;
  
  ### Now for the variances. We compute the variance in a normal
  ### population, because when we want to test the skewness and
  ### kurtosis, we do so under the assumption that they come from a
  ### normal distribution. Plus, no idea how to compute the variance
  ### otherwise :-)
  ### These formula's come from the Joanes & Gill (1998) article.
  res$intermediate$var.g1 <- ( 6 * (n-2) ) /
                             ( (n+1) * (n+3) );
  res$intermediate$var.g2 <- ( 24*n * (n-2) * (n-3) ) /
                             ( (n+1)^2 * (n+3) * (n + 5) );
  res$intermediate$var.G1 <- res$intermediate$var.g1 *
                             sqrt(n * (n-1)) /
                             (n-2)^2;
  res$intermediate$var.G2 <- res$intermediate$var.g2 *
                             (((n-1) * (n+1)) /
                              ((n-2) * (n-3))) ^ 2;
  res$intermediate$var.b1 <- res$intermediate$var.g1 *
                             ((n-1) / n) ^ 3;
  res$intermediate$var.b2 <- res$intermediate$var.g2 *
                             ((n-1) / n) ^ 4;
  
  ### And the standard error for the skewness and the kurtosis,
  ### based on the formulas on http://www.tc3.edu/instruct/sbrown/stat/shape.htm
  res$intermediate$se.G1 <- sqrt( ( 6 * n * (n-1) ) /
                                  ( (n-2) * (n+1) * (n+3) ) );
  res$intermediate$se.G2 <- 2 * res$intermediate$se.G1 *
                            sqrt( (n^2 - 1) /
                                  ( (n-3) * (n+5) ) );
  
  ### Confidence interval for skewness and kurtosis in the population
  res$intermediate$zMultiplier <- qnorm(1 - (1-conf.level)/2);
  res$intermediate$ci.G1 <- c(res$intermediate$G1 -
                                res$intermediate$zMultiplier *
                                res$intermediate$se.G1,
                              res$intermediate$G1 +
                                res$intermediate$zMultiplier *
                                res$intermediate$se.G1);
  res$intermediate$ci.G2 <- c(res$intermediate$G2 -
                                res$intermediate$zMultiplier *
                                res$intermediate$se.G2,
                              res$intermediate$G2 +
                                res$intermediate$zMultiplier *
                                res$intermediate$se.G2);
  
  ### Z-test to test against normal distribution
  res$intermediate$z.G1 <- res$intermediate$G1 / res$intermediate$se.G1;
  res$intermediate$z.G2 <- res$intermediate$G2 / res$intermediate$se.G2;
  res$intermediate$p.G1 <- 2*(1-pnorm(abs(res$intermediate$z.G1)));
  res$intermediate$p.G2 <- 2*(1-pnorm(abs(res$intermediate$z.G2)));
  
  ### Store requested estimates in output object
  if (type == 1) {
    res$output$skewness <- res$intermediate$g1;
    res$output$kurtosis <- res$intermediate$g2;
    res$output$type <- "g";
  }
  else if (type == 2) {
    res$output$skewness <- res$intermediate$G1;
    res$output$kurtosis <- res$intermediate$G2;
    res$output$type <- "G";
  }
  else if (type == 3)  {
    res$output$skewness <- res$intermediate$b1;
    res$output$kurtosis <- res$intermediate$b2;
    res$output$type <- "b";
  }
  
  ### Return results object
  class(res) <- "dataShape";
  return(res);
  
}

print.dataShape <- function(x, digits=x$input$digits, extraNotification=TRUE, ...) {
  if (x$output$type == "G") {
    skewness.inference <- paste0("  (se = ", format(x$intermediate$se.G1, digits=digits),
                                 ", confidence interval = [",
                                 paste0(format(x$intermediate$ci.G1, digits=digits), collapse=", "),
                                 "], z = ", format(x$intermediate$z.G1, digits=digits),
                                 ", p = ", format(x$intermediate$p.G1, digits=digits), ")");
    kurtosis.inference <- paste0("  (se = ", format(x$intermediate$se.G2, digits=digits),
                                 ", confidence interval = [",
                                 paste0(format(x$intermediate$ci.G2, digits=digits), collapse=", "),
                                 "], z = ", format(x$intermediate$z.G2, digits=digits),
                                 ", p = ", format(x$intermediate$p.G2, digits=digits), ")");
  }
  else {
    skewness.inference <- kurtosis.inference <- "";
  }
  cat(paste0("Skewness (", x$output$type, "1): ",
             round(x$output$skewness, digits=digits),
             skewness.inference, "\n"));
  cat(paste0("Kurtosis (", x$output$type, "2): ",
             round(x$output$kurtosis, digits=digits),
             kurtosis.inference, "\n"));
  if (extraNotification && x$output$type == "g") {
    cat("\nNote: g1 and g2 are biased estimates of the population skewness and kurtosis.",
        "For unbiased estimates, use 'type=2' or 'type=3'.");
  }
  else if (extraNotification && x$output$type == "G") {
    cat("\nNote: G1 and G2 are the estimates for skewness and kurtosis used by SPSS and SAS,",
        "and corrected for the bias present in g1 and g2 ('type=1'). Note that b1 and b2 ('type=3')",
        "may perform better in small samples from a normal distribution.");
  }
  else if (extraNotification && x$output$type == "b") {
    cat("\nNote: b1 and b2 are estimates for skewness and kurtosis that have smaller mean-squared error in small",
        "samples from a normal distribution than G1 and G2, which are used by SPSS and SAS ('type=2').");
  }
}


### Function definition
normalityAssessment <- function(sampleVector, samples = 5000, digits=3,
                                samplingDistColor = "#2222CC",
                                normalColor = "#00CC00",
                                samplingDistLineSize = 2,
                                normalLineSize = 1,
                                xLabel.sampleDist = NULL,
                                yLabel.sampleDist = NULL,
                                xLabel.samplingDist = NULL,
                                yLabel.samplingDist = NULL) {
  
  ### Create object for returning results
  res <- list(sampleVector.raw = sampleVector,
              sampleVector = sampleVector[complete.cases(sampleVector)],
              sampleSize = length(sampleVector[complete.cases(sampleVector)]),
              samples = samples,
              digits = digits);
  
  ### Construct temporary dataset for
  ### plotting sample distribution
  normalX <- c(seq(min(res$sampleVector), max(res$sampleVector),
                   by=(max(res$sampleVector) - min(res$sampleVector))/(res$sampleSize-1)));
  normalY <- dnorm(normalX, mean=mean(res$sampleVector),
                   sd=sd(res$sampleVector));
  sampleDistY <- res$sampleVector;
  tempDat <- data.frame(normalX = normalX, normalY = normalY, sampleDist = sampleDistY);
  tempBinWidth <- (max(res$sampleVector) - min(res$sampleVector)) / 30;

  ### Generate labels if these weren't specified
  if (is.null(xLabel.sampleDist)) {
    xLabel.sampleDist <- paste0('Value of ', deparse(substitute(sampleVector)));
  }
  if (is.null(yLabel.sampleDist)) {
    yLabel.sampleDist <- paste0('Frequencies for n=', res$sampleSize);
  }
  
  ### Plot sample distribution
  res$plot.sampleDist <- powerHist(tempDat$sampleDist,
                                   xLabel=xLabel.sampleDist,
                                   yLabel=yLabel.sampleDist,
                                   distributionColor=samplingDistColor,
                                   normalColor=normalColor,
                                   distributionLineSize=samplingDistLineSize,
                                   normalLineSize=normalLineSize);
  
  ### Take 'samples' samples of sampleSize people and store the means
  ### (first generate an empty vector to store the means)
  res$samplingDistribution <- c();
  for (i in 1:samples) {
    res$samplingDistribution[i] <- mean(sample(res$sampleVector, size=res$sampleSize,
                                               replace=TRUE));
  }
  
  ### Construct temporary dataset for
  ### plotting sampling distribution  
  normalX <- c(seq(min(res$samplingDistribution), max(res$samplingDistribution),
                   by=(max(res$samplingDistribution) - min(res$samplingDistribution))/(res$samples-1)));
  normalY <- dnorm(normalX, mean=mean(res$samplingDistribution),
                   sd=sd(res$samplingDistribution));
  samplingDistY <- res$samplingDistribution;
  tempDat <- data.frame(normalX = normalX, normalY = normalY, samplingDist = samplingDistY);
  tempBinWidth <- (max(res$samplingDistribution) - min(res$samplingDistribution)) / 30;
  
  ### Generate labels if these weren't specified
  if (is.null(xLabel.samplingDist)) {
    xLabel.samplingDist <- paste0('Value of ', deparse(substitute(sampleVector)));
  }
  if (is.null(yLabel.samplingDist)) {
    yLabel.samplingDist <- paste0('Frequencies for ', res$samples, ' samples of n=', res$sampleSize);
  }

  ### Plot sampling distribution
  res$plot.samplingDist <- powerHist(tempDat$samplingDist,
                                     xLabel=xLabel.samplingDist,
                                     yLabel=yLabel.samplingDist,
                                     distributionColor=samplingDistColor,
                                     normalColor=normalColor,
                                     distributionLineSize=samplingDistLineSize,
                                     normalLineSize=normalLineSize);
  
  ### Shapiro Wilk test - if there are more than 5000
  ### datapoints, only use the first 5000 datapoints
  res$sw.sampleDist <- ifelseObj(res$sampleSize > 5000,
                              shapiro.test(res$sampleVector[1:5000]),
                              shapiro.test(res$sampleVector));
  res$sw.samplingDist <- ifelseObj(res$samples > 5000,
                                shapiro.test(res$samplingDistribution[1:5000]),
                                shapiro.test(res$samplingDistribution));
  
  ### Anderson-Darling test
  res$ad.sampleDist <- adTest(res$sampleVector);
  res$ad.samplingDist <- adTest(res$samplingDistribution);
  
  ### Kolomogorov-Smirnof test
  suppressWarnings(res$ks.sampleDist <-
                     ks.test(res$sampleVector, "pnorm", alternative = "two.sided"));
  suppressWarnings(res$ks.samplingDist <-
                     ks.test(res$samplingDistribution, "pnorm", alternative = "two.sided"));

  ### Skewness and kurtosis
  res$dataShape.sampleDist <- dataShape(res$sampleVector);
  res$dataShape.samplingDist <- dataShape(res$samplingDistribution);

  ### Set class for returnable object and return it
  class(res) <- 'normalityAssessment';
  return(res);
  
}

print.normalityAssessment <- function (x, ...) {

  if (x$sampleSize > 5000) {
    sw.sampleDist <- paste0("Shapiro-Wilk: p=", round(x$sw.sampleDist$p.value, x$digits),
                                   " (W=", round(x$sw.sampleDist$statistic, x$digits),
                                   "; NOTE: based on the first 5000 of ",
                                 x$sampleSize, " observations)");
  }
  else {
    sw.sampleDist <- paste0("Shapiro-Wilk: p=", round(x$sw.sampleDist$p.value, x$digits),
                                   " (W=", round(x$sw.sampleDist$statistic, x$digits),
                                   "; based on ", x$sampleSize, " observations)");
  }
  
  if (x$samples > 5000) {
    sw.samplingDist <- paste0("Shapiro-Wilk: p=", round(x$sw.samplingDist$p.value, x$digits),
                       " (W=", round(x$sw.samplingDist$statistic, x$digits),
                       "; NOTE: based on the first 5000 of ",
                       x$samples, " observations)");
  }
  else {
    sw.samplingDist <- paste0("Shapiro-Wilk: p=", round(x$sw.samplingDist$p.value, x$digits),
                           " (W=", round(x$sw.samplingDist$statistic, x$digits),
                           "; based on ", x$samples, " observations)");
  }
  
  ### Show output
  cat("## SAMPLE DISTRIBUTION ###\n");
  cat(paste0("Sample distribution of ", x$sampleSize,
             " observations\n",
             "Mean=", round(mean(x$sampleVector), x$digits),
             ", median=", round(median(x$sampleVector), x$digits),
             ", SD=", round(sd(x$sampleVector), x$digits),
             ", and therefore SE of the mean = ",
             round(sd(x$sampleVector)/sqrt(x$sampleSize), x$digits),
             "\n\n"));
  print(x$dataShape.sampleDist, extraNotification=FALSE);
  cat(paste0("\n", sw.sampleDist, "\n",
             "Anderson-Darling: p=", round(x$ad.sampleDist@test$p.value, x$digits),
             " (A=", round(x$ad.sampleDist@test$statistic, x$digits), ")\n",
             "Kolmogorov-Smirnof: p=", round(x$ks.sampleDist$p.value, x$digits),
             " (D=", round(x$ks.sampleDist$statistic, x$digits), ")"));
  
  cat("\n\n## SAMPLING DISTRIBUTION FOR THE MEAN ###\n");
  cat(paste0("Sampling distribution of ", x$samples, " samples of n=", x$sampleSize, "\n",
             "Mean=", round(mean(x$samplingDistribution), x$digits),
             ", median=", round(median(x$samplingDistribution), x$digits),
             ", SD=", round(sqrt(var(x$samplingDistribution)), x$digits),
             "\n\n"));
  print(x$dataShape.samplingDist, extraNotification=FALSE);
  cat(paste0("\n", sw.samplingDist, "\n",
             "Anderson-Darling: p=", round(x$ad.samplingDist@test$p.value, x$digits),
             " (A=", round(x$ad.samplingDist@test$statistic, x$digits), ")\n",
             "Kolmogorov-Smirnof: p=", round(x$ks.samplingDist$p.value, x$digits),
             " (D=", round(x$ks.samplingDist$statistic, x$digits), ")"));

  ### Plots
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2)));
  suppressWarnings(print(x$plot.sampleDist$plot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1)));
  suppressWarnings(print(x$plot.samplingDist$plot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2)));

  invisible(); 
}