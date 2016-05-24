chi2CorrAge <- function(formula, data.obs, namepara1, namepara2, nameage, w1, w2, mort, a, nsimu, nbcore = 3) {
  
  ## check valididy
  if(!is.numeric(w1) | !is.numeric(w2))
    stop("'w1' and 'w2' must be numeric")
  if(nsimu <= 0 | is.na(nsimu) | !is.numeric(nsimu))
    stop("'nsimu' must be a positive number")
  if(!is.vector(mort) | !is.vector(a))
    stop("'mort' and 'a' must be numeric")
  if(length(a) != length(mort) + 1)
    stop("'mort' must have the length of 'a' plus 1")
  if(class(formula) != "character")
    stop("'formula' must be a string of character")
  if(!is.character(namepara1) | !is.character(namepara2) | !is.character(nameage))
    stop("'namepara1', 'namepara2' and 'nameage' must be strings of character")
  if(length(data.obs[[nameage]]) == 0)
    stop("'nameage' must be a column name in 'data'")
  if(!all(dim(data.obs) == dim(na.omit(data.obs)))) {
    warning("data containing NA are omitted")
    ## rows having some NA values are removed
    data.obs <- na.omit(data.obs)
  }
  
  ## Time value at the beginning of calculations
  begintime <- Sys.time()
  
  
  ## Observed data
  obs <- obsdata_chi2corrage(formula = formula, data = data.obs, name1 = namepara1, name2 = namepara2, nameage = nameage, w1 = w1, w2 = w2, mort = mort, a = a, v0para1 = NA, v0para2 = NA)
  chi2.corr.obs <- obs$chi2corrobs
  infect.proba <- obs$infectproba
  
  ## with parallel computing
  if(nbcore > 1) {
    ## create the cluster
    if(requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)) {
      cl <- parallel::makeCluster(nbcore)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, c("calcInfectProba", "SensTransMatrix", "EstimParam", "ModelClass"))  ## functions in 'internal_functions.R'
      
      ## Simulated data
      
      chi2.corr.sim <- foreach::"%dopar%"(foreach::foreach(simu = 1:nsimu, .combine = c), 
                                          simudata_chi2corrage(formula = formula, data = data.obs, name1 = namepara1, name2 = namepara2, nameage = nameage, w1 = w1, w2 = w2, mort = mort, a = a, v0para1 = infect.proba$estimparam1, v0para2 = infect.proba$estimparam2, matprobainfect = infect.proba$matprobainfect))
      chi2.corr.sim <- unlist(chi2.corr.sim)
      pval <- round((1 + sum(chi2.corr.sim > chi2.corr.obs)) / (1 + nsimu), digits = 5)
        
      ## close the cluster
      parallel::stopCluster(cl)
    }
    
  ## without parallel computing
  } else {
    ## Simulated data
    ## chi2corrsim is a vector in which there is the chi-square between the simulated data and the theoretical infectius status of each Simulated individual
    chi2.corr.sim <- numeric(length = nsimu)
    for(nn in 1:nsimu) {
      chi2.corr.sim[nn] <- simudata_chi2corrage(formula = formula, data = data.obs, name1 = namepara1, name2 = namepara2, nameage = nameage, w1 = w1, w2 = w2, mort = mort, a = a, v0para1 = infect.proba$estimparam1, v0para2 = infect.proba$estimparam2, matprobainfect = infect.proba$matprobainfect)
    }
    chi2.corr.sim <- unlist(chi2.corr.sim)
    pval <- round((1 + sum(chi2.corr.sim > chi2.corr.obs)) / (1 + nsimu), digits = 5)
  }
  
  
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
  return(list(formula = formula, time = time, nbcore = nbcore, chi2.corr.obs = chi2.corr.obs, pval = pval, tab.th = obs$tabth, tab.obs = obs$tabobs, chi2.corr.sim = chi2.corr.sim))
}

