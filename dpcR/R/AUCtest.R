AUCtest <- function(x = x, y = y, threshold = 0.05, noise_cut = 0.05, savgol = TRUE, norm = FALSE,
                    filter.q = c(0.7, 0.8)) {
  
  # Get the raw data and assign them to a data frame containing the abscissa values (e.g., 
  # time, count, ...) and the corresponding peak value (e.g., fluorescence)
  # Pre-processing:
  # Set values below a "noise_cut" value to 0
  y[y < quantile(y, noise_cut)] <- 0
  
  # Normalise data based on the quantiles
  qval <- 0.05
  if (norm) 
    y <- (y - quantile(y, qval)) / (quantile(y, 1 - qval) - quantile(y, qval))
  
  # Smooth the data by Splines or Savitzky-Golay Smoothing (default)
  if (savgol == TRUE) {
    data <- data.frame(x, sgolayfilt(y))
  } else (
    data <- cbind(x, smooth.spline(x,y)[["yin"]])
  )
  
  # find *ALL* peaks of the input data based on findpeaks of the pracma package
  supposed_peaks <-  findpeaks(data[, 2])
  
  # Try to find out where "findpeaks" detected noise (no real peak but just a small spike), 
  # "negative peaks (no amplification, just primer dimer...)" and "positive peaks (true 
  # amplification)" based in the peak height.
  # filter.q is the quantile value of the noise and the quantile value of the negative peaks
  peak_quantiles <- c(min(supposed_peaks), 
                      quantile(supposed_peaks[, 1], filter.q), 
                      max(supposed_peaks))
  
  # test_res is checked later which element of the data is TRUE
  test_res <- cut(supposed_peaks[, 1], peak_quantiles, include.lowest = TRUE)
  levels(test_res) <- c("noise", "negative", "positive")
  
  sp <- smooth.spline(x, y)
  # create an empty matrix with the results of the area under the curve calculation
  # the column number of the matrix might grow depending addition of further elements
  all_peaks <- data.frame(do.call("rbind", lapply(1L:nrow(supposed_peaks), function(i) {
    # select range of single peak
    xy <- data[supposed_peaks[i, 3]:supposed_peaks[i, 4], 1:2]
    
    peak <- rep(0, 6)
    
    psp <- function(x, sp) 
      predict(sp, x)[["y"]]
    
    # Estimate the AUC by integration. NOTE: Needs improvements because the integration will 
    # fail badly if peak overlap considerably!
    try(integrate_tmp <- integrate(psp, lower = min(xy[, 1]), upper = max(xy[, 1]), sp)$value)
    peak[1] <- supposed_peaks[i, 2] # Position of the peak maximum
    #could use quadinf{pracma} or adaptIntegrate{cubature} instead of integrate, need further investigation
    peak[2] <- ifelse(class(integrate_tmp) == "try-error", NA, integrate_tmp) # Crude estimation of the AUC
    peak[3] <- max(xy[, 1]) - min(xy[, 1]) # crude estimation of the peak width
    peak[4] <- supposed_peaks[i, 1] # height of the peak depending on time ...
    peak[5] <- data[supposed_peaks[i,2],1] # position of the peak depending on the time ...
    # Determine the resolution of the data
    peak[6] <- mean(vapply(2L:nrow(xy), function(i) abs(xy[i, 1] - xy[i - 1, 1]), 0))  # time resolution of the peaks
    peak
  })))
  
  all_peaks <- cbind(1L:nrow(supposed_peaks), test_res, all_peaks)
  colnames(all_peaks) <- c("Peak number", "State", "Position", "AUC", "Width", "Height", 
                           "Index", "Resolution")
  
  list(peaks = all_peaks, data = data) # List containing the table with all values
  # and the smoothed data
}