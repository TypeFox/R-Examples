## -----------------------------------------------------------------------------
## Molecular diffusion coefficients
## -----------------------------------------------------------------------------

## Boudreau 1997 - Arrhenius formulation
## H2, CH4, He, Ne, Kr, Xe, Rn: Jahne et al (1987)
## DMS: Saltzman et al. (1993)
## Ar: Ohsumi and Horibe (1984)

.marelac$ArrDat <- data.frame(
  A  = c(3338, 3047, 2000, 818, 1608, 7238, 6393, 9007, 15877),
  Ea = c(16.06,18.36,18.10,11.70,14.84,19.81,20.20,21.61,23.26)
)
rownames(.marelac$ArrDat) <- c("H2","CH4","DMS","He","Ne","Ar","Kr","Xe","Rn")

## from Wilke and Chang (1955) as modified by Hayduk and Laudie (1974)
.marelac$WCDat <- data.frame(A=c(34.7,35.2,24.5,23.6,36.0,34.5,43.8),
                             B = NA)
rownames(.marelac$WCDat) <- c("N2","H2S","NH3","NO","N2O","CO","SO2")

## Boudreau (1997): linear functions of temp
.marelac$BDat <- data.frame(
  m0 = c(c(25.9,10.0,9.60,6.29,9.81,5.06,4.33,4.02,3.26,2.62,10.4,6.35,
            4.82,5.99,4.88,4.66,10.3,9.50,54.4,4.43,6.06,6.06,10.3,7.82,
            9.50,3.60,3.43,3.31,3.18,4.06,2.57,3.31,3.31,3.39,3.63,3.36,
            3.69,4.46,3.91,3.31,2.79,2.95,2.78,2.71)),
  m1 = c(1.094,0.441,0.438,0.343,0.432,0.275,0.199,0.223,
         0.177,0.143,0.273,0.280,0.266,0.307,0.232,0.252,
         0.331,0.388,1.555,0.241,0.297,0.297,0.416,0.359,
         0.413,0.179,0.144,0.150,0.155,0.176,0.140,0.152,
         0.152,0.158,0.208,0.130,0.169,0.198,0.199,0.151,
         0.172,0.131,0.136,0.120))

rownames(.marelac$BDat) <- c("OH","Br","Cl","F","I","HCO3","CO3","H2PO4",
    "HPO4","PO4","HS","HSO3","SO3","HSO4","SO4","IO3","NO2","NO3",
    "H","Li","Na","K","Cs","Ag","NH4","Ca","Mg","Fe","Mn","Ba","Be",
    "Cd","Co","Cu","Hg","Ni","Sr","Pb","Ra","Zn","Al","Ce","La","Pu")


## -----------------------------------------------------------------------------
## Molecular diffusion coefficients
## -----------------------------------------------------------------------------

diffcoeff <- function(S = 35, t = 25, P = 1.013253, species = c("H2O",
  "O2", "CO2", "H2", "CH4", "DMS", "He", "Ne", "Ar", "Kr", "Xe", "Rn",
  "N2", "H2S", "NH3", "NO", "N2O", "CO", "SO2", "OH", "F", "Cl", "Br",
  "I", "HCO3", "CO3", "H2PO4", "HPO4", "PO4", "HS", "HSO3", "SO3",
  "HSO4", "SO4", "IO3", "NO2", "NO3", "H", "Li", "Na", "K", "Cs",
  "Ag", "NH4", "Ca", "Mg", "Fe", "Mn", "Ba", "Be", "Cd", "Co", "Cu",
  "Hg", "Ni", "Sr", "Pb", "Ra", "Zn", "Al", "Ce", "La", "Pu", "H3PO4",
  "BOH3", "BOH4", "H4SiO4")) {

  ## thpe: check length of arguments to avoid ambigous behavior
  len <- c(length(S), length(t), length(P))
  if (! all(len %in% c(1, max(len))))
    stop("vectors of S, t and P do not match")

  if (any (S < 0)) stop ("Salinity should be >= 0")
  
  species <- match.arg(species, several.ok = TRUE)
  diffc   <- list()                         # will have the result
  diffArr <- diffChang <- diffBoud <- NULL  # initialize variables
  TK      <- t + 273.15                     # temperature, Kelvin
  Patm    <- 1.013253

  ##  The viscosity in pure water at atmospheric pressure and sample
  ##  temperature.
  mu_0 <- viscosity(S = 0, t = t, P = Patm)

  ##  Diffusion coefficient of Water at S = 0, P = P=1.013253, t = t
  ##  Cohen MH and Turnbull D (1959).
  ##  Krynicki K, Green CD and Sawyer DW (1978).
  if ("H2O" %in% species) {
    A         <- 12.5e-09 * exp(-5.22e-04 * P)
    B         <- 925.0 * exp(-2.6e-04 * P)
    T0        <- 95.0 + 2.61e-02 * P
    D_H2O     <- A * sqrt(TK) * exp(-B / (TK - T0))
    diffc$H2O <- D_H2O * viscosity(S = 0, t, P) / mu_0
  }

  ##  Diffusion coefficient of O2 and CO2
  ##  Boudreau (1997)
  if ("O2" %in% species) {
    A <- 0.2604
    B <- 0.006383
    diffc$O2 <- (A + B * (TK / mu_0)) * 1E-09
  }

  if ("CO2" %in% species) {
    A <- 0.1954
    B <- 0.005089
    diffc$CO2 <- (A + B * (TK / mu_0)) * 1E-09
  }

  ## H3PO4 : Least (1984) determined D(H3PO4) at 25 deg C and 0 S.
  ##         Assume that this value can be scaled by the Stokes-Einstein
  ##         relationship to any other temperature.
  if ("H3PO4" %in% species) {
    D_H3PO4 <- 0.87e-09
    tS      <- 25.0
    SS      <- 0.0
    mu_S    <- viscosity(SS, tS, Patm)
    diffc$H3PO4 <- D_H3PO4 * (mu_S / mu_0) * (TK / (tS + 273.15))
  }

  #  B(OH)3 : Mackin (1986) determined D(B(OH)3) at 25 deg C and
  #           about 29.2 S.
  #           Assume that this value can be scaled by the Stokes-Einstein
  #           relationship to any other temperature.
  if ("BOH3" %in% species) {
    D_BOH3 <- 1.12e-09
    tS     <- 25.0
    SS     <- 29.2
    mu_S   <- viscosity(SS, tS, Patm)
    diffc$BOH3 <- D_BOH3 * (mu_S / mu_0) * (TK / (tS + 273.15))
  }

  ##  B(OH)4 : No information on this species ! Boudreau and
  ##           Canfield (1988) assume it is 12.5% smaller than B(OH)3.
  if ("BOH4" %in% species) {
    if (is.null(diffc$BOH3)) {
      D_BOH3 <- 1.12e-09
      tS     <- 25.0
      SS     <- 29.2
      mu_S   <- viscosity(SS, tS, Patm)
      D_BOH3 <- D_BOH3 * (mu_S / mu_0) * (TK / (tS + 273.15))
    } else D_BOH3 <- diffc$BOH3
    diffc$BOH4 <- 0.875 * D_BOH3
  }

  ##  H4SiO4 : Wollast and Garrels (1971) found D(H4SiO4) at 25 deg C
  ##           and 36.1 ppt S.
  ##           Assume that this value can be scaled by the Stokes-Einstein
  ##           relationship to any other temperature.
  if ("H4SiO4" %in% species) {
    D_H4SiO4 <- 1.0E-09
    tS    <- 25.0
    SS    <- 36.1
    mu_S  <- viscosity(SS, tS, Patm)
    diffc$H4SiO4 <- D_H4SiO4 * (mu_S / mu_0) * (TK / (tS + 273.15))
  }

  ##  Other dissolved substances
  ##  Boudreau (1997)
  ##  TK <- 298.15
  Arrhenius <- function(A) A[1] * exp(-(A[2] * 1000)/(8.314472 * TK)) * 1.0E-09

  ii  <- .marelac$ArrDat[species[which(species %in% rownames(.marelac$ArrDat))],]
  if (nrow(ii) > 0) diffArr <- apply(ii, 1, Arrhenius)

  ##  Other dissolved substances
  ##  from Wilke and Chang (1955) as modified by Hayduk and Laudie (1974)
  WilkeChang <- function(Vb) 4.72E-07 * TK / (mu_0 * Vb[1]^0.6) * 1.0E-04

  ii <- .marelac$WCDat[species[which(species %in% rownames(.marelac$WCDat))],]
  if (length(ii) > 0) diffChang <- apply(data.frame(ii), 1, WilkeChang)

  ##  The coefficients in pure water for the following species are
  ##  calculated by linear functions of temperature (deg C)
  ##  coefficients as in Boudreau (1997).
  Boudreau <- function(m) (m[1] + m[2] * t) * 1.0e-10
  ii <- .marelac$BDat[species[which(species %in% rownames(.marelac$BDat))],]
  if (nrow(ii) > 0) diffBoud <- apply(ii, 1, Boudreau)

  ## Simplify this
  if (is.matrix(diffc)    | is.matrix(diffChang) |
      is.matrix(diffBoud) | is.matrix(diffArr)) {
    Diffc <- NULL
           
    ## thpe: xcbind is a non-exported function with more argument checking
    Diffc <- xcbind(Diffc, diffc)
    Diffc <- xcbind(Diffc, diffArr)
    Diffc <- xcbind(Diffc, diffChang)
    Diffc <- xcbind(Diffc, diffBoud)
    
    
    ## thpe: alternative function with double check
    ##       remove this after testing phase
    #if (length(diffc)     > 0) Diffc <- xcbind(Diffc, as.data.frame(diffc))
    #if (length(diffArr)   > 0) Diffc <- xcbind(Diffc, diffArr)
    #if (length(diffChang) > 0) Diffc <- xcbind(Diffc, data.frame(diffChang))
    #if (length(diffBoud)  > 0) Diffc <- xcbind(Diffc, diffBoud)
    
    Diffc <- as.data.frame(Diffc)
  } else Diffc <- data.frame(c(diffc, diffArr, diffChang, diffBoud))

  ## SALINITY AND PRESSURE CORRECTION
  ##  To correct for pressure and salinity, the Stokes-Einstein relationship
  ##  is used. This is not quite accurate, but it is at least consistent.

  mu <- viscosity(S, t, P)   #  viscosity at sample conditions
  fac <- (mu_0 / mu)

  ## Ks -> ThPe HERE AGAIN, if salinity is a vector, then fac is also a vector
  ## and I use a loop here !
  if (nrow(Diffc) != length(fac)) {
    DD <- Diffc
    for (i in 1:(length(fac) - 1)) DD <- rbind(DD, Diffc)
    Diffc <- DD
  }
  diffc <- Diffc * fac

  return(diffc[species])  # [species] to have same ordering as input
}
