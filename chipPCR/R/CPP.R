CPP <- function(x, y, smoother = TRUE, method = "savgol", trans = FALSE, 
                method.reg = "lmrob", bg.outliers = FALSE, median = FALSE, 
                method.norm = "none", qnL = 0.03, amptest = FALSE, manual = FALSE, 
                nl = 0.08, bg.range = NULL, ...) {
  
  testxy(x, y)
  
  # Define the method for the smoothing
  # functionality identical to smoother function

  
  method.smooth <- check.method(c("savgol", "lowess", "mova", "smooth", "spline", 
                 "supsmu", "whit1", "whit2"), method)
  
  # Remove missing values from y
  if(any(is.na(y)))
    y <- fixNA(x, y, spline = TRUE)
  
  #Try to detect background range automatically or manually
  if (is.null(bg.range)){
    bg <- bg.max(x, y)
    BG <- slot(bg, "bg.start"):slot(bg, "bg.stop")
  } else {
    if(length(bg.range) != 2)
      stop("bg.range must have two values - start and end of the background.")
    BG <- min(bg.range):max(bg.range)
  }
  # Remove outliers (high and low) from the background range and 
  # substitute with the median
  if (bg.outliers) { 
    y[BG] 	<- rm.outlier(y[BG], fill = TRUE, 
                         median = median)
    y[BG] 	<- rm.outlier(y[BG], opposite = TRUE, 
                         fill = TRUE, median = median)
  }
  
  # Invoke smoother to improve data for further calculations
  # SERVE BUG (Priority high): "mova" will caus problems because it 
  # truncates the data (mova 3 -> first and last value miss at the end)
  
  if (smoother) {
    y <- smoother(x, y, trans = FALSE, bg.outliers = FALSE, 
                  method = list(method.smooth), CPP = FALSE, ...)
    #some smoothing methods (for example mova) can introduce missing values
    if(any(is.na(y))) {
      y <- fixNA(x, y, spline = TRUE)
      message(paste0("NA values introduced by smoothing method ", 
                     method.smooth, ". fixNA() used."))
    }
  } else {
    y
  }
  
  # Test if linear correction based on the background range is requested
  # If requested first try a robust linear regression. If robust linear 
  # regression fails try standard lm()
  if (trans) {
    coefficients <- lm.coefs(x[BG], y[BG], method.reg)
    # Apply linear model to the raw data
    y.norm <- y - (coefficients[2, 1] * x + coefficients[1, 1])
    # Subtract the median (based on background range) from the data
    y.norm <- y.norm - median(y.norm[BG]) 
    # Subtract the median (based on background range) from the data 
    # without a linear model
  } else {y.norm <- y - median(y[BG])}
  
  # Perform a normalization to a specified quantile value
  y.norm <- normalizer(y = y.norm, method.norm = method.norm, qnL = qnL)
 
  # Test if the amplifification is likely to be positive
  if (amptest) {
    y.norm <- amptester(y.norm, manual = manual, 
                        background = range(BG), 
                        noiselevel = nl)
  }
  
  list(y.norm = y.norm, BG = BG)
}


setGeneric("CPP")


setMethod("CPP", signature(x = "data.frame", y="missing"), 
          function(x, y, smoother = TRUE, trans = FALSE, method.reg = "lmrob", 
                   bg.outliers = FALSE, median = FALSE, 
                   method.norm = "none", qnL = 0.03, 
                   amptest = FALSE, manual = FALSE, nl = 0.08, bg.range = NULL, ...) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            CPP(x[, 1], x[, 2], trans = trans, method.reg = "lmrob", ,
                bg.outliers = bg.outliers, 
                median = median, method.norm = method.norm, qnL = qnL,
                amptest = amptest, manual = manual, nl = nl, bg.range = NULL, ...)
          })

setMethod("CPP", signature(x = "matrix", y="missing"), 
          function(x, y, smoother = TRUE, trans = FALSE, method.reg = "lmrob", 
                   bg.outliers = FALSE, median = FALSE, 
                   method.norm = "none", qnL = 0.03, 
                   amptest = FALSE, manual = FALSE, nl = 0.08, bg.range = NULL, ...) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            CPP(x[, 1], x[, 2], smoother = TRUE, trans = trans, 
                method.reg = "lmrob", bg.outliers = bg.outliers, 
                median = median, method.norm = method.norm, qnL = qnL,
                amptest = amptest, manual = manual, nl = nl, bg.range = NULL, ...)
          })
