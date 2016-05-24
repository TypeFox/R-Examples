## -----------------------------------------------------------------------------
## Salinity-Chlorinity Conversion
## -----------------------------------------------------------------------------

convert_StoCl <- function(S=35) {
  if (any (S<0))
    stop ("Salinity should be >= 0")

  S / 1.80655 # chlorinity, in g/kg
}
