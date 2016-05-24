# Modified 5 Nov 2015 SDH

vol2mol <- function(
  vol,           # Volume in mL
  gas = 'CH4',
  temp,
  pres,
  rh = 0,
  unit.temp = getOption('unit.temp', 'C'),
  unit.pres = getOption('unit.pres', 'atm'),
  tp.message = TRUE
) {

  gas <- toupper(gas)

  checkArgClassValue(gas, expected.values = names(vol.mol))

  if(tp.message) message('You specified a gas pressure of ', pres, ' ', unit.pres, ' and temperature of ', temp, ' ', unit.temp, '.')

  # Convert to standard volume (101.325 kPa, 273.15 K, dry)
  pres.pa <- unitConvert(x = pres, unit = unit.pres, to = 'Pa')
  temp.k <- unitConvert(x = temp, unit = unit.temp, to = 'K')
  vol <- stdVol(vol, temp = temp.k, pres = pres.pa, rh = rh, pres.std = 101325, temp.std = 273.15, unit.temp = 'K', unit.pres = 'Pa', std.message = FALSE)

  # Gas quantity in moles from standard volume
  mol <- vol/vol.mol[gas]

  return(mol)

}
