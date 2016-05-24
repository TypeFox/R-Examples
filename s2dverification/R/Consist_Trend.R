Consist_Trend <- function(var_exp, var_obs, interval = 1) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 6 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 3 | length(dimobs) < 3) {
    stop("At least 3 dim needed : c(nexp/nobs, nsdates, nltime)") 
  }
  for (jn in 2:max(length(dimexp), length(dimobs))) {  
    if (dimexp[jn] != dimobs[jn]) {
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 6)
  var_obs <- Enlarge(var_obs, 6)
  #
  #  Find common points to compute trends 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  tmp <- MeanListDim(var_obs, dims = 4:6, narm = TRUE)
  tmp2 <- MeanListDim(var_exp, dims = 4:6, narm = TRUE)
  nan <- Mean1Dim(tmp, 1, narm = FALSE) + Mean1Dim(tmp2, 1, narm = FALSE)
  for (jdate in 1:dimexp[2]) {
    for (jtime in 1:dimexp[3]) {
      if (is.na(nan[jdate, jtime])) {
        var_exp[, jdate, jtime, , , ] <- NA
        var_obs[, jdate, jtime, , , ] <- NA
      }
    }
  }
  #
  #  Compute trends 
  # ~~~~~~~~~~~~~~~~
  #
  trendcat <- array(dim = c(dimobs[1] + dimexp[1], 4, dim(var_obs)[3:6]))
  detrendedmod <- array(dim = dimexp)
  detrendedobs <- array(dim = dimobs)
  tmp1 <- Trend(var_exp, 2, interval = interval)
  tmp2 <- Trend(var_obs, 2, interval = interval)
  trendcat[1:dimexp[1], , , , , ] <- tmp1$trend
  trendcat[(dimexp[1] + 1):dim(trendcat)[1], , , , , ] <- tmp2$trend
  detrendedmod[] <- tmp1$detrended
  detrendedobs[] <- tmp2$detrended
  #
  #  Back to the original dimensions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  outtrendcat <- array(dim = c(dimobs[1] + dimexp[1], 4, 
                               dimobs[3:length(dimobs)]))
  outtrendcat[] <- trendcat
  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(trend = outtrendcat, detrendedmod = detrendedmod,
                 detrendedobs = detrendedobs))
}
