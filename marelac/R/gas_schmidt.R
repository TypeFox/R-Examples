## -----------------------------------------------------------------------------
## The Schmidt Number for a Gas
## -----------------------------------------------------------------------------

## Pre-load local data
.marelac$SchmidtCoeff <- data.frame(
  A = c(410.14,855.1,2206.1,1638,1909.1,2205,3412.8,2039.2,2073.1,
        2301.1,3845.4,3501.8,3531.6,4295.8),
  B = c(20.503,46.299,144.86,81.83,125.09,135.71,224.3,120.31,125.62,
        151.1,228.95,210.31,231.4,281.52),
  C = c(0.53175,1.254,4.5413,1.483,3.9012,3.9549,6.7954,3.4209,3.6276,
        4.7364,6.1908,6.1851,7.2168,8.7826),
  D = c(0.006011,0.01449,0.056988,0.008004,0.048953,0.047339,0.083,0.040437,
        0.043219,0.059431,0.06743,0.07513,0.090558,0.11025)
)

rownames(.marelac$SchmidtCoeff) <- c("He", "Ne", "N2", "O2", "Ar",
        "Kr", "Rn", "CH4","CO2", "N2O", "CCl2F2", "CCL3F",
        "SF6", "CCl4")


## The function itself
## thpe: make the original version to an internal one 'gas_schmidt1'
gas_schmidt1 <- function (t = 25, species = c("He", "Ne", "N2", "O2", "Ar",
        "Kr", "Rn", "CH4","CO2", "N2O", "CCl2F2", "CCL3F", "SF6", "CCl4")) {

  species <- match.arg(species, several.ok = TRUE)
  Sc <- .marelac$SchmidtCoeff[species,]
  schmidt <- Sc$A - Sc$B * t + Sc$C * t*t - Sc$D * t*t*t
  schmidt <- t(matrix (nrow = length(species), data = schmidt))
  colnames(schmidt) <- species
  schmidt
}

## ThPe: now define the user-visible vectorized version
gas_schmidt <- function (t = 25, species = c("He", "Ne", "N2", "O2", "Ar",
        "Kr", "Rn", "CH4","CO2", "N2O", "CCl2F2", "CCL3F", "SF6", "CCl4")) {
        
  if ((length(t) == 1)) {
    ## return vector for species OR temperatures
    ret <- gas_schmidt1(t, species)
  } else {
    ## return a matrix for temperatures x species
    ret <- t(sapply(X = t, FUN = function(X) gas_schmidt1(t = X, species)))
    if (is.matrix(ret) & min(dim(ret)) > 1) 
      colnames(ret) <- species
    else
      ret <- as.vector(ret)
  }
  ret
}

