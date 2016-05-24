amptester <-
  function(y, manual = FALSE, noiselevel = 0.08, background = NULL) {
    testxy(x = y, y, both = FALSE)
    # Test if background has only two values
    if (!is.null(background) && length(background) != 2)
      stop("Use only two values (e.g., background = c(1,10)) \n\t to set the range for the background correction")
    if (is.null(background) && manual == TRUE)
      stop("Manual test requires specified background.")
    #if background is NULL, sorting it is pointless and invokes warning
    if (!is.null(background))
      background <- as.integer(sort(background))
    
    # fix possible missing vaues with fixNA (spline method)
    y <- fixNA(1L:length(y), y)
    
    # FIRST TEST
    # Shapiro test (SHt)
    # Simple test if data come from noise or presumably a melting curve
    noisy <- FALSE
    res.shapiro <- shapiro.test(y)[["p.value"]]
    if (res.shapiro >= 5e-04) {
      message("The distribution of the curve data indicates noise.\n\t The data should be visually inspected.")
      noisy <- TRUE
      mess.shapiro <- "Appears not to be an amplification curve"
    } else {
      mess.shapiro <- "Appears to be an amplification curve"
    }
    
    # Determine the head and tail region for the next tests.
    # Apply a simple rule to take the first 20 percent and the last 15 percent
    # of any input data set to calculate the number of elements for the head 
    # (nh) and tail (nt), to deal with other data types like isothermal
    # amplifications
    
    nh <- trunc(length(y) * 0.2)
    if (nh < 5) nh <- 5
    
    nt <- trunc(length(y) * 0.15)
    if (nt < 5) nt <- 5
    
    # SECOND TEST
    # Resids growth test (RGt)
    # test if fluorescence values in linear phase are stable. Whenever no amplification 
    # occurs, fluorescence values quickly deviate from linear model. Their standarized
    # residuals will be strongly correlated with their value. For real amplification curves,
    # situation is much more stable. Noise (that means deviations from linear model) 
    # in  background do not correlate strongly with the changes in fluorescence. 
    # The decision is based on the threshold value (here 0.5). 
    rgts <-sapply(0L:round(length(y)/8, 0), function (j) {
      cyc <- 1:nh + j
      reg <- lm(y[cyc] ~ cyc)
      cor(rstudent(reg), y[cyc])  
    })
    rgt.dec <- ifelse(sum(rgts < 0.8) > round(length(y)/16, 0), "positive", "negative")
    
    # THIRD TEST
    # Linear Regression test (LRt)
    # This test determines the R^2 by a linear regression. The R^2 are
    # determined from a run of circa 15 percent range of the data.
    # If a sequence of more than six R^2s is larger than 0.8 is found 
    # thant is likely a nonlinear signal. This is a bit counterintuitive 
    # because R^2 of nonlinear data should be low.
    
    ws <- ceiling((15 * length(y)) / 100)
    if (ws < 5) 
      ws <- 5
    if (ws > 15) 
      ws <- 15
    y.tmp <- na.omit(y[-c(1:5)])
    x <- 1:length(y.tmp)
    suppressWarnings(
      res.reg <- sapply(1L:(length(y.tmp)), function (i)  {
        round(summary(lm(y.tmp[i:c(i + ws)] ~ x[i:c(i + ws)]))[["r.squared"]], 4)
      }
      )
    )
    
    # Binarize R^2 values. Everything larger than 0.8 is positve
    res.LRt <- res.reg
    # Define the limits for the R^2 test
    res.LRt[res.LRt < 0.8] <- 0
    res.LRt[res.LRt >= 0.8] <- 1
    # Seek for a sequence of at least six positve values (R^2 >= 0.8)
    # The first five measurepoitns of the amplification curve are skipped
    # beacuse most technologies and probetechnologies tend to overshot
    # in the start (background) region.
    res.out <- sapply(5L:(length(res.LRt) - 6), function(i) {
      ifelse(sum(res.LRt[i:(i + 4)]) == 5, TRUE, FALSE)
    })
    
    # Test if more than one sequence of positive values was found (will 
    # be the case in most situation due to an overlap of the positive 
    # sequences.)
    linearity <- ifelse(sum(res.out) >= 1, TRUE, FALSE)
    
    # FOURTH TEST (MANUAL)
    # Threshold test (THt)
    # Manual test for positve amplification based on a fixed threshold
    # value.
    if (manual) {
      signal <- median(y[-(background)]) - mad(y[-(background)])
      if (signal < noiselevel && (signal / noiselevel) < 1.25) {
        y <- abs(rnorm(length(y), 0, 0.1^30))
        tht.dec <- "negative"
      } else {
        tht.dec <- "positive"
      }
    } else {
      # FOURTH TEST (AUTOMATIC)
      # Threshold test (THt)
      # Apply a simple rule to take the first 20 percent and the last 15 percent
      # of any input data set and perform a Wilcoxon rank sum tests for the head 
      # (nh) and tail (nt).
      
      res.wt <- suppressWarnings(wilcox.test(head(y, n = nh), tail(y, n = nt), 
                                             alternative = "less"))
      
      if (res.wt$p.value > 0.01) {
        y <- abs(rnorm(length(y), 0, 0.1^30))
        tht.dec <- "negative"
      } else {
        tht.dec <- "positive"
      }
    }
    
    # FIFTH TEST
    # Signal level test (SLt)
    # The meaninfulness can be tested by comparison of the signals
    # 1) A robust "sigma" rule by median + 2 * mad 
    # 2) comparison of the signal/noise ratio. If less than 1.25 (25 percent) 
    # signal increase it is likely that nothing happened during the reaction.
    noisebackground <- median(head(y, n = nh)) + 2 * mad(head(y, n = nh))
    signal  <- median(tail(y, n = nt)) - 2 * mad(tail(y, n = nt))
    if (signal <= noisebackground || signal / noisebackground <= 1.25) {
      y <- abs(rnorm(length(y), 0, 0.1^30))
      slt.dec <- "negative"
    } else {
      slt.dec <- "positive"
    }
    
    # SIXTH TEST
    # The pco test determines if the points in an amplification curve (like a polygon)
    # are in a "clockwise" order. The sum over the edges result in a positive value if the
    # amplification curve is "clockwise" and is negative if the curve is counter-clockwise.
    # From experience is noise positive and "true" amplification curves "highly" negative.
    # This test depends on the definition of a threshold.
    pco <- function(y) {
      xy <- data.frame(predict(smooth.spline(1L:length(y), y)))
      sum(sapply(1L:(nrow(xy) - 1), 
                 function (i) {
                   xy[i + 1, 1] - xy[i, 1] * xy[i + 1, 2] + xy[i, 2]
                 })
      )
    }
    
    der.res <- summary(inder(1L:length(y), y), print = FALSE)
    
    lm.dat <- data.frame(x = c(round(der.res[["SDM"]], 0), round(der.res[["SDm"]], 0)))
    lm.dat <- cbind(lm.dat, y = y[lm.dat[, 1]])      
    lm.dat[["y"]] <- lm.dat[["y"]]/max(lm.dat[["y"]])
    slope.ratio <- coef(lm(y ~ x, lm.dat))
    res.pco <- pco(y)
    
    # Output of the different tests
    rgt.dec <- ifelse(rgt.dec == "positive", TRUE, FALSE)
    tht.dec <- ifelse(tht.dec == "positive", TRUE, FALSE)
    slt.dec <- ifelse(slt.dec == "positive", TRUE, FALSE)
    new("amptest", 
        .Data = y, 
        decisions = c(shap.noisy = noisy,
                      lrt.test = linearity,
                      rgt.dec = rgt.dec,
                      tht.dec = tht.dec,
                      slt.dec = slt.dec), 
        noiselevel = noiselevel,
        background = background,
        polygon = res.pco,
        slope.ratio = slope.ratio[["x"]])
  }
