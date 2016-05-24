# Modified: 11 Nov 2015 SDH

vol2mass <- function(
  volBg,
  xCH4,
  temp.hs,
  temp.vol,
  pres.hs,
  pres.vol,
  unit.temp = getOption('unit.temp', 'C'),
  unit.pres = getOption('unit.pres', 'atm'),
  rh.hs = 1,
  rh.vol = 1
) {

  # Check arguments
  checkArgClassValue(volBg, c('numeric', 'integer'), expected.range = c(0, Inf))
  checkArgClassValue(xCH4, c('numeric', 'integer'), expected.range = c(0, 1))
  checkArgClassValue(temp.hs, c('numeric', 'integer'))
  checkArgClassValue(temp.vol, c('numeric', 'integer'))
  checkArgClassValue(pres.hs, c('numeric', 'integer'))
  checkArgClassValue(pres.vol, c('numeric', 'integer'))
  checkArgClassValue(unit.temp, c('character'))
  checkArgClassValue(unit.pres, c('character'))
  checkArgClassValue(rh.hs, c('numeric', 'integer'), expected.range = c(0, 1))
  checkArgClassValue(rh.vol, c('numeric', 'integer'), expected.range = c(0, 1))

  # Standardize volBg
  # Must be standardized to 0C and 101325 Pa, because molar volumes are for these conditions
  # First convert temperature and pressures to K and Pa
  temp.hs.k <- unitConvert(x = temp.hs, unit = unit.temp, to = 'K')
  pres.hs.pa <- unitConvert(x = pres.hs, unit = unit.pres, to = 'Pa')
  temp.vol.k <- unitConvert(x = temp.vol, unit = unit.temp, to = 'K')
  pres.vol.pa <- unitConvert(x = pres.vol, unit = unit.pres, to = 'Pa')

  volBg <- stdVol(volBg, temp = temp.vol.k, pres = pres.vol.pa, rh = rh.vol, temp.std = 273.15, pres.std = 101325, unit.temp = 'K', unit.pres = 'Pa', std.message = TRUE)

  # Calculate molar mass and molar volume
  # Volumes defined at 101325 Pa (1 atm) and 273.15 K (0C)
  mmb <- xCH4*molMass('CH4') + (1 - xCH4)*molMass('CO2')
  mvb <- xCH4*vol.mol['CH4'] + (1 - xCH4)*vol.mol['CO2']
  mvb <- as.vector(mvb)
  # Dry, standardized biogas density
  db <- mmb/mvb

  # Calculate saturated water vapor pressure in atm (based on NIST)
  pH2O <- rh.hs*watVap(temp.k = temp.hs.k) 

  # And mass of water in biogas
  mH2O <- molMass('H2O')*pH2O/((pres.hs.pa - pH2O)*mvb)

  # Mass loss in g
  mass <- volBg*(db + mH2O)

  return(mass)
}
