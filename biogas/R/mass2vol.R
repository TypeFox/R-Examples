# Modified: 20 Feb 2016 SDH

mass2vol <- function(
  mass,
  xCH4,
  temp,
  pres,
  temp.std = getOption('temp.std', 0.0),
  pres.std = getOption('pres.std', 1.0),
  unit.temp = getOption('unit.temp', 'C'),
  unit.pres = getOption('unit.pres', 'atm'),
  value = "CH4",
  headspace = NULL,
  headcomp = 'N2',
  temp.init = NULL, # = unitConvert(20, unit = 'C', to = unit.temp),
  std.message = TRUE
) {

  # Hardwire rh, maybe will add as argument again later
  rh <- 1

  # Check arguments~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  checkArgClassValue(mass, c('integer', 'numeric'), expected.range = c(0, Inf))
  checkArgClassValue(xCH4, c('integer', 'numeric'), expected.range = c(0, 1))
  checkArgClassValue(temp, c('integer', 'numeric'))
  checkArgClassValue(pres, c('integer', 'numeric'))
  checkArgClassValue(temp.std, c('integer', 'numeric'))
  checkArgClassValue(pres.std, c('integer', 'numeric'))
  checkArgClassValue(unit.temp, c('character'))
  checkArgClassValue(unit.pres, c('character'))
  checkArgClassValue(value, c('character'))
  checkArgClassValue(tolower(value), expected.values = c('ch4', 'co2', 'bg', 'all'))
  checkArgClassValue(headspace, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(headcomp, c('character', 'NULL'))
  checkArgClassValue(temp.init, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(std.message, 'logical')

  # Unit conversions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Convert pressure to atm and temp to C
  pres.pa <- unitConvert(x = pres, unit = unit.pres, to = 'Pa')
  temp.k <- unitConvert(x = temp, unit = unit.temp, to = 'K')
  if(!is.null(temp.init)) temp.init.k <- unitConvert(x = temp.init, unit = unit.temp, to = 'K')

  # Convert standard values to K and Pa
  temp.std.k <- unitConvert(x = temp.std, unit = unit.temp, to = 'K')
  pres.std.pa <- unitConvert(x = pres.std, unit = unit.pres, to = 'Pa')

  # Density calcuation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mmb <- xCH4*molMass('CH4') + (1 - xCH4)*molMass('CO2')
  mvBg <- xCH4*vol.mol['CH4'] + (1 - xCH4)*vol.mol['CO2']
  mvBg <- as.vector(mvBg)
  db <- mmb/mvBg

  # Calculate water vapor pressure in Pa (based on NIST)
  pH2O <- rh*watVap(temp.k = temp.k)

  mH2O <- molMass('H2O')*pH2O/((pres.pa - pH2O)*mvBg)

  # Biogas volume
  vBg <- mass/(db + mH2O)

  # Initial headspace correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # If headspace data provided
  if(!is.null(headspace)) {
    if(headcomp == 'N2') {
      if(is.null(temp.init)) stop('Trying to apply correction for initial headspace composition, but no value was provided for temp.init.')
      message('In gravimetric method, correcting for initial headspace composition, using an initial temperature of ', temp.init, ' ', unit.temp, '.')
      # Reduce mass loss for N2 loss
      # Fixed temperature and pressure at start and fixed residual pressure (1 atm)
      # 0.0012504 is the density of N2 in g/ml at 1C and 1 atm (101325 Pa) from NIST
      mass <- mass - stdVol(headspace, temp = temp.init.k, pres = 101325, rh = 0, temp.std = 273.15, pres.std = 101325, unit.pres = 'Pa', unit.temp = 'K', std.message = FALSE)*0.0012504
      vBg <- mass/(db + mH2O) + stdVol(headspace, temp = temp.k, pres = 101325, rh = 1, temp.std = 273.15, pres.std = 101325, unit.pres = 'Pa', unit.temp = 'K', std.message = FALSE)
    } else {
      warning('headcomp argument was given as ', headcomp, ', but \"N2\" is the only option available, so no correction was applied for initial headspace composition.')
    }
  }

  # Standardize (based on molar volumes used above, so in the default case stdVol does nothing.)
  #vBg <- stdVol(vBg, temp = unitConvert(x = 273.15, unit = 'K', to = unit.temp), pres = unitConvert(x = 101325, unit = 'Pa', to = unit.pres), rh = 0, temp.std = temp.std, pres.std = pres.std, unit.pres = unit.pres, unit.temp = unit.temp, std.message = std.message)
  vBg <- stdVol(vBg, temp = 273.15, pres = 101325, rh = 0, temp.std = temp.std.k, pres.std = pres.std.pa, unit.pres = 'Pa', unit.temp = 'K', std.message = std.message)
  vCH4 <- xCH4*vBg*vol.mol['CH4']/mvBg
  vCH4 <- as.vector(vCH4)
  vCO2 <- (1 - xCH4)*vBg*vol.mol['CO2']/mvBg
  vCO2 <- as.vector(vCO2)

  # Return output
  if(tolower(value) == "ch4") {
    return(vCH4) 
  } else if(tolower(value) == "bg") {
    return(vBg) 
  } else if(tolower(value) == "co2") {
    return(vCO2) 
  } else if(tolower(value) == "all") {
    return(cbind(vBg = vBg, vCH4 = vCH4, vCO2 = vCO2))
  }
}
