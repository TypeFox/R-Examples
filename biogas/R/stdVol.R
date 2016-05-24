# Modified: 2 April 2016 SDH

stdVol <- function(
  vol,
  temp,
  pres,
  rh = 1,
  temp.std = getOption('temp.std', 0.0),
  pres.std = getOption('pres.std', 1.0),
  unit.temp = getOption('unit.temp', 'C'),
  unit.pres = getOption('unit.pres', 'atm'),
  std.message = TRUE
) {

  # Check arguments
  checkArgClassValue(vol, c('numeric', 'integer'))
  checkArgClassValue(temp, c('numeric', 'integer'), expected.range = c(0, Inf))
  checkArgClassValue(pres, c('numeric', 'integer'), expected.range = c(0, Inf))
  checkArgClassValue(rh, c('numeric', 'integer'), expected.range = c(0, 1))
  checkArgClassValue(temp.std, c('numeric', 'integer'))
  checkArgClassValue(pres.std, c('numeric', 'integer'), expected.range = c(0, Inf))
  checkArgClassValue(unit.temp, c('character'))
  checkArgClassValue(unit.pres, c('character'))
  checkArgClassValue(std.message, c('logical'))

  # Echo standard conditions
  if(std.message) message('Using a standard pressure of ', pres.std, ' ', unit.pres, ' and standard temperature of ', temp.std, ' ', unit.temp, ' for standardizing volume.')

  # Convert pressure to Pa and temp to K
  pres.pa <- unitConvert(x = pres, unit = unit.pres, to = 'Pa')
  temp.k <- unitConvert(x = temp, unit = unit.temp, to = 'K')
  pres.std.pa <- unitConvert(x = pres.std, unit = unit.pres, to = 'Pa')
  temp.std.k <- unitConvert(x = temp.std, unit = unit.temp, to = 'K')

  # Check pres and temp range
  if(any(temp.k < 273.15 | temp.k > 373.15)) warning('temp ranges from ', min(na.omit(temp)), ' ', unit.temp, ' to ', max(na.omit(temp)), ' ', unit.temp, '. Is this really correct?')
  if(any(pres.pa < 5E4 | pres.pa > 1.5E5)) warning('pres ranges from ', min(na.omit(pres)), ' ', unit.pres, ' to ', max(na.omit(pres)), ' ', unit.pres, '. Is this really correct?')
 
  # Calculate water vapor pressure in atm (based on NIST)
  pH2O <- rh*watVap(temp.k = temp.k)

  # Correct volume for water and to standard pressure
  vol.dry <- vol*(pres.pa - pH2O)/pres.std.pa

  # Correct dry volume to standard temperature
  vol.std <- vol.dry*temp.std.k/temp.k

  return(vol.std)
}
