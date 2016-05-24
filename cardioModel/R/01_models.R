###################################################################################
######################Baseline models##############################################
###################################################################################

# Intercept (E0)
.my.intercept.fun <- function (Int) {
  Int
}

# Cosine period 24 h
.my.cosine24.fun <- function (Mean, Amp1, Tpeak1, TOD) {
  Mean+Amp1*cos(2*pi*TOD/24-Tpeak1)
}
#The optimization algorithm used in nlme does not
#allow constraints on the coefficients in the model.
#eliminate Amp and peak because etas in them are having computational difficulties

# Cosine period 12 h
.my.cosine12.fun <- function (Mean, Amp2, Tpeak2, TOD) {
  Mean+Amp2*cos(2*pi*TOD/12-Tpeak2)
}

# Double Cosine period 24 and 12 h
.my.dcosine.fun <- function (Mean, Amp1, Tpeak1, Amp2, Tpeak2, TOD) {
  Mean+Amp1*cos(2*pi*TOD/24-Tpeak1)+Amp2*cos(2*pi*TOD/12-Tpeak2)
}

# Function: Minimum time after meal (Mtam) as covariate of Intercept (E0)
.my.intercept.meal.fun <- function (Int, Mtam, Kn, Kp) {
  Int*4/((1+exp(Kn*Mtam))*(1+exp(Kp*Mtam)))
}

###################################################################################
########################Drug models################################################
###################################################################################

##Slope (linear)
.my.slope.fun <- function (Slo, Conc) {
  Slo*Conc
}

##Emax
.my.Emax.fun <- function (Emax, EC50, Conc) {
  Emax*Conc/(exp(EC50)+Conc)
}

##Sigmoidal Emax
.my.sigEmax.fun <- function (Emax, EC50, Hill, Conc) {
  Emax*Conc^(exp(Hill))/(exp(EC50)^(exp(Hill))+Conc^(exp(Hill)))
}
