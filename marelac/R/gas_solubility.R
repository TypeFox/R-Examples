## -----------------------------------------------------------------------------
## Solubility of a Gas in Seawater
## -----------------------------------------------------------------------------

## Pre-load local data

#Coefficients for the fit of solubility to the following equations:
#ln(x)=A1+A2(100/T)+A3ln(T/100)+A4(T/100)^2+S[B1+B2(T/100)+B3(T/100)^2]
#T is temperature in kelvin, S is salinity,
#In data.frame BunsenSolubCoeff (x) is the bunsen
#solubility coefficient in units of l gal (l solution)-1 atm-1,
#in data.frame volumeSolubCoeff, (x) is the volumetric solubility function,
#in units of mol /l/atm

.marelac$BunsenSolubCoeff <- data.frame(
  A1 = c(-34.6261,-39.1971,-59.6274,-58.3877,-55.6578,-57.2596,-11.95,-68.8862),
  A2 = c(43.0285,51.8013,85.7761,85.8079,82.0262,87.4242,31.66,101.4956),
  A3 = c(14.1391,15.7699,24.3696,23.8439,22.5929,22.9332,0,28.7314),
  A4 = c(0,0,0,0,0,0,0,0),
  B1 = c(-0.04234,-0.124695,-0.05158,-0.034892,-0.036267,-0.008723,0,-0.076146),
  B2 = c(0.022624,0.078374,0.026329,0.015568,0.016241,-0.002793,0,0.04397),
  B3 = c(-0.003312,-0.0127972,-0.0037252,-0.0019387,-0.0020114,
       0.0012398,0,-0.0068672),
  type = rep(1, 8),
  Accur = c(0.005,0.005,0.004,0.004,0.004,0.004,0,0.01)
)
rownames(.marelac$BunsenSolubCoeff) <-
  c("He", "Ne", "N2", "O2", "Ar", "Kr", "Rn", "CH4")

.marelac$VolumeSolubCoeff <- data.frame(
  A1 = c(-160.7333,-165.8806,-218.0971,-229.9261,-80.0343,-148.247),
  A2 = c(215.4152,222.8743,298.9702,319.6552,117.232,227.758),
  A3 = c(89.8920,92.0792,113.8049,119.4471,29.5817,62.5557),
  A4 = c(-1.47759,-1.48425,-1.39165,-1.39165,0,0),
  B1 = c(0.029941,-0.056235,-0.143566,-0.142382,0.0335183,-0.400847),
  B2 = c(-0.027455,0.031619,0.091015,0.091459,-0.0373942,0.265218),
  B3 = c(0.0053407,-0.0048472,-0.0153924,-0.0157274,0.00774862,-0.0446424),
  type = rep(2, 6),
  Accur = c(0.003,0.0014,0.015,0.015,0.02,0.025)
)
rownames(.marelac$VolumeSolubCoeff) <-
  c("CO2", "N2O", "CCl2F2", "CCl3F", "SF6", "CCl4")
.marelac$SolubCoeff <- rbind(.marelac$BunsenSolubCoeff, .marelac$VolumeSolubCoeff)

## Todo: remove redundant coefficients from .marelac

## -----------------------------------------------------------------------------
## The function itself
## -----------------------------------------------------------------------------
gas_solubility <- function (S = 35, t = 25, species = c("He", "Ne", "N2", "O2",
  "Ar", "Kr", "Rn", "CH4", "CO2", "N2O", "CCl2F2", "CCl3F", "SF6", "CCl4")) {

  if (! checkVecLength(list(S, t)))
    warning("Arguments 'S' and 't' should have the same length or length 1.")    

  if (any (S < 0)) stop ("Salinity should be >= 0")
  K <-  t + 273.15
  Sbc <- .marelac$SolubCoeff[species,]
  SA <- NULL

  for (i in 1:nrow(Sbc)) {
    Sb <- Sbc[i,]
    bet <- Sb$A1 + Sb$A2 * (100/K) + Sb$A3 * log(K/100) + Sb$A4 * (K/100)^2 +
           S * (Sb$B1 + Sb$B2 * K/100 + Sb$B3 * (K/100)^2)

    if (Sb$type == 1) SS  <- exp(bet)/22.4136 * 10^6/1.013253
    if (Sb$type == 2) SS  <- exp(bet)/1.013253 / (1 - vapor(t = t, S = S)) * 10^6
    SA <- cbind(SA, SS)
  }

  colnames(SA) <- species
  ## convert matrix to a vector if one of the dimensions is "1"
  if (min(dim(SA)) == 1) SA <- as.vector(SA)
  SA
}
