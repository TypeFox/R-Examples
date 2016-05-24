# From Hartree
Hartree2electronVolt <- function(x)
  Joule2electronVolt(Hartree2Joule(x))

Hartree2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(Hartree2Joule(x))

Hartree2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(Hartree2Joule(x))

Hartree2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(Hartree2Joule(x))

Hartree2Kelvin <- function(x)
  Joule2Kelvin(Hartree2Joule(x))

Hartree2Joule <- function(x)
  x*universalConstants["Eh","Value"]

Hartree2Hertz <- function(x)
  Joule2Hertz(Hartree2Joule(x))

# From electronVolt
electronVolt2Hartree <- function(x)
  Joule2Hartree(electronVolt2Joule(x))

electronVolt2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(electronVolt2Joule(x))

electronVolt2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(electronVolt2Joule(x))

electronVolt2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(electronVolt2Joule(x))

electronVolt2Kelvin <- function(x)
  Joule2Kelvin(electronVolt2Joule(x))

electronVolt2Joule <- function(x)
  x*universalConstants["e","Value"]

electronVolt2Hertz <- function(x)
  Joule2Hertz(electronVolt2Joule(x))

# From centimeterMinusOne
centimeterMinusOne2Hartree <- function(x)
  Joule2Hartree(centimeterMinusOne2Joule(x))

centimeterMinusOne2electronVolt <- function(x)
  Joule2electronVolt(centimeterMinusOne2Joule(x))

centimeterMinusOne2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(centimeterMinusOne2Joule(x))

centimeterMinusOne2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(centimeterMinusOne2Joule(x))

centimeterMinusOne2Kelvin <- function(x)
  Joule2Kelvin(centimeterMinusOne2Joule(x))

centimeterMinusOne2Joule <- function(x)
  x*1E2*universalConstants["h","Value"]*universalConstants["c","Value"]

centimeterMinusOne2Hertz <- function(x)
  Joule2Hertz(centimeterMinusOne2Joule(x))

# From kiloCaloriePerMole
kiloCaloriePerMole2Hartree <- function(x)
  Joule2Hartree(kiloCaloriePerMole2Joule(x))

kiloCaloriePerMole2electronVolt <- function(x)
  Joule2electronVolt(kiloCaloriePerMole2Joule(x))

kiloCaloriePerMole2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(kiloCaloriePerMole2Joule(x))

kiloCaloriePerMole2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(kiloCaloriePerMole2Joule(x))

kiloCaloriePerMole2Kelvin <- function(x)
  Joule2Kelvin(kiloCaloriePerMole2Joule(x))

kiloCaloriePerMole2Joule <- function(x)
  x*universalConstants["cal","Value"]*1E3/universalConstants["Na","Value"]

kiloCaloriePerMole2Hertz <- function(x)
  Joule2Hertz(kiloCaloriePerMole2Joule(x))

# From kiloJoulePerMole
kiloJoulePerMole2Hartree <- function(x)
  Joule2Hartree(kiloJoulePerMole2Joule(x))

kiloJoulePerMole2electronVolt <- function(x)
  Joule2electronVolt(kiloJoulePerMole2Joule(x))

kiloJoulePerMole2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(kiloJoulePerMole2Joule(x))

kiloJoulePerMole2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(kiloJoulePerMole2Joule(x))

kiloJoulePerMole2Kelvin <- function(x)
  Joule2Kelvin(kiloJoulePerMole2Joule(x))

kiloJoulePerMole2Joule <- function(x)
  x*1E3/universalConstants["Na","Value"]

kiloJoulePerMole2Hertz <- function(x)
  Joule2Hertz(kiloJoulePerMole2Joule(x))

# From Kelvin
Kelvin2Hartree <- function(x)
  Joule2Hartree(Kelvin2Joule(x))

Kelvin2electronVolt <- function(x)
  Joule2electronVolt(Kelvin2Joule(x))

Kelvin2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(Kelvin2Joule(x))

Kelvin2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(Kelvin2Joule(x))

Kelvin2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(Kelvin2Joule(x))

Kelvin2Joule <- function(x)
  x*universalConstants["kb","Value"]

Kelvin2Hertz <- function(x)
  Joule2Hertz(Kelvin2Joule(x))

# From Joule
Joule2Hartree <- function(x)
  x/universalConstants["Eh","Value"]

Joule2electronVolt <- function(x)
  x/universalConstants["e","Value"]

Joule2centimeterMinusOne <- function(x)
  x*1E-2/(universalConstants["h","Value"]*universalConstants["c","Value"])

Joule2kiloCaloriePerMole <- function(x)
  Joule2kiloJoulePerMole(x)/universalConstants["cal","Value"]

Joule2kiloJoulePerMole <- function(x)
  x*universalConstants["Na","Value"]*1E-3

Joule2Kelvin <- function(x)
  x/universalConstants["kb","Value"]

Joule2Hertz <- function(x)
  x/universalConstants["h","Value"]

# From Hertz
Hertz2Hartree <- function(x)
  Joule2Hartree(Hertz2Joule(x))

Hertz2electronVolt <- function(x)
  Joule2electronVolt(Hertz2Joule(x))

Hertz2centimeterMinusOne <- function(x)
  Joule2centimeterMinusOne(Hertz2Joule(x))

Hertz2kiloCaloriePerMole <- function(x)
  Joule2kiloCaloriePerMole(Hertz2Joule(x))

Hertz2kiloJoulePerMole <- function(x)
  Joule2kiloJoulePerMole(Hertz2Joule(x))

Hertz2Kelvin <- function(x)
  Joule2Kelvin(Hertz2Joule(x))

Hertz2Joule <- function(x)
  x*universalConstants["h","Value"]
