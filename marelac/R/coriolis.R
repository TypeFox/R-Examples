## -----------------------------------------------------------------------------
## Coriolis Factor as a Function of Latitude
## -----------------------------------------------------------------------------

coriolis <- function (lat)  # latitude in degrees north (-90:+90)
  2 * 7.2921e-5 * sin(pi/180 * lat)
