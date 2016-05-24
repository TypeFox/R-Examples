MFIaggr <- function(x, y, cyc = 1, fluo = 2:ncol(x), RSD = FALSE, 
		    rob = FALSE, llul = c(1,10)){
  
  #cyc and fluo arguments are used in methods
  #Test if x and y exist.
  testxy(x, y, length = FALSE)
  
  # Test if y has enough values
  if (is.null(dim(y)))
    warning("MFIaggr does not work properly with small number of fluorescence data columns.")
  
  
  # Test if llul has only two values
  if (!is.null(llul) && length(llul) != 2)
    stop("Use two cycle values (e.g., llul = c(1,10)) to set 
        the range for the analysis.")
  
  # llul defines the range of interrest (ROI) to be analyzed
  # in detail. In particular, early cycles (background range)
  # and final cycles (plateau phase) are potential targets.
  llul <- sort(llul)	
  
  # Decide which method for the calculation of location and 
  # dispersion should be used. "rob" means robust and it will
  # use median and MAD instead of mean and sd
  
  if (rob) {
    loc.fct <- median
    dev.fct <- mad
  } else {
    loc.fct <- mean
    dev.fct <- sd
  }
  
  # NOTE: no error hanling in if only on collumn of input data is
  # used. Waring: apply will fail in this case. Need fix: yes, see
  # effcalc
  
  # Use apply over the input data to calculate the location and 
  # dispersion of a data bulk
  
  y.m <- apply(y, 1, loc.fct)
  y.sd <- apply(y, 1, dev.fct)
  
  # Decide if the relative standard deviation should be 
  # calculated
  
  if (RSD) {
    y.cv <- (y.sd / y.m) * 100
  } else {
    y.cv <- y.sd / y.m
  }
  
  # Apply the results to the data.frame "res" and label
  # collumns according to the used location and dispersion
  # method
  res <- data.frame(x, y.m, y.sd, y.cv)
  
  if (rob) {
    if(RSD) {
      names(res)  <- c("Cycle", "Location (Median)", 
                       "Deviation (MAD)", 
                       "Coefficient of Variance (RSD)")
    } else {
      names(res) <- c("Cycle", "Location (Median)", 
                      "Deviation (MAD)", 
                      "Coefficient of Variance (RSD [%])")
    }
  } else {
    if(RSD) {
      names(res)  <- c("Cycle", "Location (Mean)", 
                       "Deviation (SD)", 
                       "Coefficient of Variance (RSD)")
    } else {
      names(res)  <- c("Cycle", "Location (Mean)", 
                       "Deviation (SD)", 
                       "Coefficient of Variance (RSD [%])")
    }
  }
  
  # Trend of the ROI
  # Calculate the trend of the ROI by a linear function. The
  # input values for the y values are calculated from the
  # mean or median depending on the setting of rob.
  
  fluo <- res[c(llul[1]:llul[2]), 2]
  cycles <- res[c(llul[1]:llul[2]), 1]
  
  lm.roi <- lm(fluo ~ cycles)
  summ.lm <- summary(lm.roi)

  # Calcuate robust und non-robust location and dispersion
  # parameters of the ROI and apply the results to stats
  y.roi <- na.omit(unlist(y[llul[1]:llul[2], ]))
  
  mean.roi <- mean(y.roi)
  median.roi <- median(y.roi)
  sd.roi <- sd(y.roi)
  
  # test for heteroscedasticity
  heter.p <- bptest(lm.roi)[["p.value"]][["BP"]]
  
  stats <- c(mean = mean.roi,
	     median = median.roi, 
	     sd = sd.roi, 
	     mad = mad(y.roi), 
	     IQR = IQR(y.roi), 
	     medcouple = mc(y.roi), 
	     skewness = 3* (mean.roi - median.roi) / sd.roi,
	     SNR = mean.roi / sd.roi, 
	     VRM = var(y.roi) / mean.roi,
	     NAs = sum(is.na(y.roi)),
	     intercept = summ.lm[["coefficients"]][1],
	     slope = summ.lm[["coefficients"]][2],
	     r.squared = summ.lm[["r.squared"]],
	     heter.p = heter.p
  )
  
  res.dens <- density(y.roi)
  res.qq <- qqnorm(y.roi, plot.it = FALSE)
  #res is the an object of the type data.frame containing the 
  #temperature, location, deviation and coefficient of variance.
  new("refMFI", .Data = res, density = res.dens, 
      qqnorm.data = as.data.frame(y[llul, ]), stats = stats, summ.lm = summ.lm)  
}

setGeneric("MFIaggr")


setMethod("MFIaggr", signature(x = "data.frame", y="missing"), 
          function(x, y, cyc = 1, fluo = 2:ncol(x), RSD = FALSE, 
		   rob = FALSE, llul = c(1,10)) { 
            MFIaggr(x[, cyc], x[, fluo], RSD = RSD, rob = rob, llul = llul)
          })

setMethod("MFIaggr", signature(x = "matrix", y="missing"), 
          function(x, y, cyc = 1, fluo = 2:ncol(x), RSD = FALSE, 
		   rob = FALSE, llul = c(1,10)) { 
            MFIaggr(x[, cyc], x[, fluo], RSD = RSD, rob = rob, llul = llul)
          })

