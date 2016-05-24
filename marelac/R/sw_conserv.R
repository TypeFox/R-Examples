## -----------------------------------------------------------------------------
## concentrations of conservative species in seawater
## -----------------------------------------------------------------------------

sw_conserv<- function(S = 35) {
  if (any (S < 0))
    stop ("Salinity should be >= 0")
  Borate   <- 4.16e2 * S / 35
  Calcite  <- 0.01028e6 * S / 35
  Sulphate <- convert_StoCl(S) * 0.14e6 * 1 / 96.062
  Fluoride <- convert_StoCl(S) * 67 / 18.9984
  return(data.frame(Borate   = Borate,
                    Calcite  = Calcite,
                    Sulphate = Sulphate,
                    Fluoride = Fluoride))
}
