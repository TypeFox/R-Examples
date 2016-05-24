## -----------------------------------------------------------------------------
## Atmospheric Composition
## -----------------------------------------------------------------------------

## Pre-load local data
.marelac$atmComp <- c(He = 5.24e-06, Ne = 1.818e-05, N2 = 0.78084, O2 = 0.20946,
    Ar = 0.00934, Kr = 1.14e-06, CH4 = 1.745e-06, CO2 = 0.000365,
    N2O = 3.14e-07, H2 = 5.5e-07, Xe = 8.7e-08, CO = 5e-08, O3 = 1e-08)

## The function
atmComp <- function(species = c("He", "Ne", "N2", "O2", "Ar", "Kr", "CH4",
                        "CO2", "N2O", "H2", "Xe", "CO", "O3")) {
  species <- match.arg(species, several.ok = TRUE)
  .marelac$atmComp[species]
}

