## -----------------------------------------------------------------------------
## Convert from Mol to Liter For a Gas
## -----------------------------------------------------------------------------

## Pre-load local data
.marelac$Waals.a <- c(0,1.363, 3.640, 11.77, 1.505, 19.7483, 6.579, 13.04, 12.18,
      20.19, 4.692, 2.283, 9.649, 19.26, 8.779, 5.536, 0.03457,
      0.2476, 4.510, 3.716, 4.490, 8.200, 2.349, 4.225, 0.2135,
      1.358, 1.408, 5.354, 3.832, 1.378, 4.692, 4.377, 4.251, 6.803, 4.250)

.marelac$Waals.b <- c(0, 0.03219, 0.04267, 0.07685, 0.03985, 0.1281,  0.05622, 0.09213,
      0.08407, 0.1286,  0.05264, 0.04278, 0.06702, 0.146,   0.08445, 0.03049,
      0.0237, 0.02661, 0.04431, 0.04081, 0.04287, 0.01696, 0.03978, 0.03707,
      0.01709, 0.02789, 0.03913, 0.04424, 0.04415, 0.03183, 0.05156, 0.05786,
      0.05571, 0.05636, 0.05105)

## The function itself
molvol <- function(t=25,
                  P=1.013253,
                  species=c("ideal","Ar", "CO2", "CS2", "CO", "CCl4", "Cl2",
      "C2H6S", "C2H5OH","C6H5F", "CH3F", "CH4", "CH3OH", "C5H12", "C3H8",
      "H2O", "He","H2", "HBr", "HCl", "H2S", "Hg", "Kr", "NH3", "Ne",
      "NO", "N2", "NO2", "N2O", "O2", "PH3", "SiH4", "SiF4", "SO2", "Xe"),
                  quantity=1, a=0,  b=0) {
  Names <- eval(formals(sys.function(sys.parent()))$species)
  ## The following is needed if called by other function, e.g. "outer";
  ## is there a better way to do this avoiding doubled information?
  Names <- c("ideal","Ar", "CO2", "CS2", "CO", "CCl4", "Cl2",
      "C2H6S", "C2H5OH","C6H5F", "CH3F", "CH4", "CH3OH", "C5H12", "C3H8",
      "H2O", "He","H2", "HBr", "HCl", "H2S", "Hg", "Kr", "NH3", "Ne",
      "NO", "N2", "NO2", "N2O", "O2", "PH3", "SiH4", "SiF4", "SO2", "Xe")
  if (! is.null(species)) {
    species <- match.arg(species, several.ok = TRUE) # check if valid input...
    ii <- pmatch(species,Names) # position of species
    # a in L2bar/mol2
    aa <- .marelac$Waals.a[ii]
    bb <- .marelac$Waals.b[ii]
  } else {
    aa<-a
    bb<-b
  }
#  R  <- 0.082058*1.0131253    #  l*bar/K/mol
  R  <- 0.0831447215   #  l*bar/K/mol

  TK <- 273.15 + t  #  t in degrees Kelvin

  il <- max(length(t), length(P), length(quantity))
  if (il > 1) {
    TK       <- rep(TK, len = il)
    P        <- rep(P, len = il)
    quantity <- rep(quantity, len = il)
  }

  VV <- NULL
  what <-  NULL
  for (i in 1:length(aa)) {
    a <- aa[i]
    b <- bb[i]

    V <- NULL
    for (i in 1:il) {
     TT <- TK[i]
     PP <- P[i]
     xx <- quantity[i]
          if (a==0 & b==0) {
            V    <- c(V, xx * R * TT / PP)   # CHECK FACTOR 1.013253
          } else {
            V <- c(V,
            uniroot(function (V) ((PP + xx * xx * a/(V^2)) *
              (V/xx - b) - R * TT), c(-10, 1e6))$root)
      }
    }
    VV <- cbind(VV, V)
  }
  colnames(VV) <- species
  if(nrow(VV) == 1) {VV <- as.vector(VV); names(VV) <- species}
  return(VV)
}

