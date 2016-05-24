## -----------------------------------------------------------------------------
## Saturated Concentration of Gasses in Seawater
## -----------------------------------------------------------------------------

gas_satconc <- function(S = 35, t = 25, P = 1.013253,
          species = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4", "CO2", "N2O"),
          atm = atmComp(species)) {
  if (any (S < 0))
    stop ("Salinity should be >= 0")
  if (! checkVecLength(list(S, t, P)))
    warning("Arguments 'S', 't' and 'P' should have the same length or length 1.")    

  if (is.null(atm)) atm <- atmComp(species)
  Vapor <- vapor(t = t, S = S)
  if (length(Vapor) == 1)
    gas_solubility(S = S, t = t, species = species) * P * atm * (1 - Vapor)
  else {
    gs <- gas_solubility(S = S, t = t, species = species)
    if (is.matrix(gs)) {
      gs * P * matrix(data = atm, nrow = nrow(gs), ncol = ncol(gs), byrow = TRUE) *
        matrix(data = (1 - Vapor), nrow = nrow(gs), ncol = ncol(gs))
    } else {
      gas_solubility(S = S, t = t, species = species) * P * atm * (1 - Vapor)    
    }
  }
}
