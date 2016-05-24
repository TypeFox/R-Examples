# Modified: 222 JUNE 2015 SDH
# NTS: would be good to switch to Pa and K for base units

unitConvert <- function(x, unit, to) {

  if(unit==to) return(x)

  # Pressure, first convert to atm
  # Temperature, first convert to C
  if(unit=="kPa") {
  x <- x/101.325 
  } else if(unit=="hPa") {
  x <- x/1013.25
  } else if(unit=="Pa") {
  x <- x/1.01325E5
  } else if(unit=="bar") {
  x <- x/1.01325
  } else if(unit=="psi") {
  x <- x/14.69595  
  } else if(unit=="atm") {
  x <- x
  } else if(unit=="F") {
  x <- (x - 32)*5/9 
  } else if(unit=="K") {
  x <- x - 273.15  
  } else if(unit=="C") {
  x <- x
  } else {
    stop('\"unit\" argument not recognized. \nOptions are \"atm\", \"kPa\", \"hPa\", \"Pa\", \"bar\", \"psi\" for pressure and \n\"C\", \"F\", or \"K\" for temperature')
  }

  # Then convert to "to" units
  if(to=="kPa") {
    x <- x*101.325 
  } else if(to=="hPa") {
    x <- x*1013.25
  } else if(to=="Pa") {
    x <- x*1.01325E5
  } else if(to=="bar") {
    x <- x*1.01325
  } else if(to=="psi") {
    x <- x*14.69595  
  } else if(to=="atm") {
  x <- x
  } else if(to=="F") {
  x <- x*9/5 + 32
  } else if(to=="K") {
  x <- x + 273.15  
  } else if(to=="C") {
  x <- x
  } else {
    stop('\"to\" argument not recognized. \nOptions are \"atm\", \"kPa\", \"hPa\", \"Pa\", \"bar\", \"psi\" for pressure and \n\"C\", \"F\", or \"K\" for temperature')
  }

  # Check coherence of conversion (at the end so that an unidentified unit will return the correct error above). 
  p.units = c('kPa', 'hPa', 'Pa', 'bar', 'psi', 'atm')
  t.units = c('F', 'K', 'C')
  if(!(unit %in% p.units & to %in% p.units) & !(unit %in% t.units & to %in% t.units)) stop('Conversion mixes temperature and pressure units.')

  return(x) 
}
