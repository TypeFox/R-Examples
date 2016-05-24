chi2Corr <- function(formula, data.obs, namepara1, namepara2, nsimu) {
   
  ## check valididy
  if(nsimu <= 0 | is.na(nsimu) | !is.numeric(nsimu))
    stop("'nsimu' must be a positive number")
  if(class(formula) != "character")
    stop("'formula' must be a string of character")
  if(!is.character(namepara1) | !is.character(namepara2))
    stop("'namepara1' and 'namepara2' must be strings of character")
  if(!all(dim(data.obs) == dim(na.omit(data.obs)))) {
    warning("data containing NA are omitted")
    ## rows having some NA values are removed
    data.obs <- na.omit(data.obs)
  }

  ## Time value at the beginning of calculations
  begintime <- Sys.time()
  
  ## Observed data 
  obs <- obsdata_chi2corr(formula = formula, data = data.obs, name1 = namepara1, name2 = namepara2)
  chi2.corr.obs <- obs$chi2corrobs
  
  ## Simulated data
  simu <- simudata_chi2corr(formula = formula, data = data.obs, name1 = namepara1, name2 = namepara2, nbsimu = nsimu, pvir1 = obs$pvir1, pvir2 = obs$pvir2, chi2corrobs = chi2.corr.obs)
  chi2.corr.sim <- simu$chi2corrsim
  
  ## Time value at the end of calculations
  endtime <- Sys.time()
  
  ## Calculation time
  time <- round(difftime(endtime, begintime, units = "secs"), digits = 3)
  
  ## Plot
 	if(chi2.corr.obs > max(chi2.corr.sim))
    maxx <- ceiling(chi2.corr.obs)
  else
    maxx <- ceiling(max(chi2.corr.sim))
  
  hist(chi2.corr.sim, xlim = c(0, maxx), main = "Position of the observed value in the distribution of the bootstrapped corrected Chi2", xlab = "Chi2corr", tcl = 0.35, cex.main = 0.9, cex = 0.7)
 	box(which = "plot", lty = "solid")
 	points(chi2.corr.obs, 0, type = "p", pch = "*", col = "red3", bg = "red3", cex = 4)
 
  ## Print in the R consol
  return(list(formula = formula, time = time, chi2.corr.obs = chi2.corr.obs, dispcoeff = simu$dispcoeff, pval1 = simu$pval1, pval2 = simu$pval2, tab.th = obs$tabth, tab.obs = obs$tabobs, chi2.corr.sim = chi2.corr.sim))
}

