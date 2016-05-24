BchronCalibrate <-
function(ages,ageSds,calCurves,ids=NULL,positions=NULL,pathToCalCurves=system.file('data',package='Bchron'),eps=1e-5,dfs=rep(100,length(ages))) {
  
  # This function expects ages in years BP (either 14C or not depending on calCurve values)
  # and positions (usually depths) in cm
  # It scales these by ageScaleVal and positionScaleVal respectively to calibrate and then
  # returns them in the original units
  
  # Check lengths of everything
  if(length(ages)!=length(ageSds)) stop("ages and ageSds should be of same length")
  if(length(ages)!=length(calCurves)) stop("ages and calCurves should be of same length")
  if(!is.null(positions)) if(length(ages)!=length(positions)) stop("ages and positions should be of same length")
  if(is.null(ids)) ids = paste('date',1:length(ages),sep='')
  
  # Check that ages and ageSds are whole numbers (i.e. years)
  if(!all(as.integer(ages)==ages)) {
    ages = round(ages,0)
    warning("ages not given as whole numbers - rounding occurred")
  }
  if(!all(as.integer(ageSds)==ageSds)) {
    # Smallest sd is 1
    ageSds = pmax(round(ageSds,0),1)
    warning("ageSds not given as whole numbers - rounding occurred")
  }
  
  # Load in all calibration curves specified
  allCalCurves = unique(calCurves) 
  calCurve = calBP = c14BP = calSd = ageGrid = mu = tau1 = list()
  for(i in 1:length(allCalCurves)) {
    calCurveFile = paste(pathToCalCurves,'/',allCalCurves[i],'.txt.gz',sep='')
    if(!file.exists(calCurveFile)) stop(paste('Calibration curve file',calCurveFile,'not found'))
    calCurve = as.matrix(utils::read.table(calCurveFile))
    calBP[[i]] = calCurve[,1]
    c14BP[[i]] = calCurve[,2]
    calSd[[i]] = calCurve[,3]
    # Create an age grid and get mean and variance of calibration curve 
    ageGrid[[i]] = seq(min(calBP[[i]]),max(calBP[[i]]),by=1)
    mu[[i]] = stats::approx(calBP[[i]],c14BP[[i]],xout=ageGrid[[i]])$y
    tau1[[i]] = stats::approx(calBP[[i]],calSd[[i]],xout=ageGrid[[i]])$y
    # Allow for greater age ranges if the calibration curve is normal
    if(allCalCurves[i]=='normal') {
      ageRange = range(c(calBP,ages+4*ageSds))
      ageGrid[[i]] = seq(ageRange[1],ageRange[2],by=1)
      mu[[i]] = ageGrid[[i]]
      tau1[[i]] = rep(0,length(ageGrid[[i]]))
    }
    
  }
  matchCalCurves = match(calCurves,allCalCurves)
      
  # Storage
  out = list()
  
  # Loop through ages and calculate densities
  for(i in 1:length(ages)) {
    
    # Get rid of ages outside the range of the calibration curve
    if(ages[i]>max(ageGrid[[matchCalCurves[i]]]) | ages[i]<min(ageGrid[[matchCalCurves[i]]])) stop(paste("Date",ids[i],"outside of calibration range"))
    
    tau = ageSds[i]^2 + tau1[[matchCalCurves[i]]]
    
    currAgeGrid = ageGrid[[matchCalCurves[i]]]
    dens = stats::dt((ages[i]-mu[[matchCalCurves[i]]])/sqrt(tau),df=dfs[i])
    dens = dens/sum(dens)
  
    # Create list of output
    if(is.null(positions)) {
      out[[i]] = list(ages=ages[i],ageSds=ageSds[i],calCurves=calCurves[i],ageGrid=currAgeGrid[dens>eps],densities=dens[dens>eps])
    } else {
      out[[i]] = list(ages=ages[i],ageSds=ageSds[i],positions=positions[i],calCurves=calCurves[i],ageGrid=currAgeGrid[dens>eps],densities=dens[dens>eps])
      
    }
  }

  names(out) = ids  
  class(out) = 'BchronCalibratedDates'
  return(out)
  
}
