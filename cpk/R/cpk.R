# Oscar A. Linares MD and David T. Daly JD/MBA
# Plymouth Pharmacokinetic Modeling Study Group
# Plymouth, Michigan 48170
#
# cpk Version 1.3-1
# 12/24/2013

#####################################
# R FUNCTIONS
#####################################

# target therapeutic concentration
ttc.fn <- function(msc, mec) {
  numerator   <- msc - mec
  denominator <- log(msc/mec)
  result      <- numerator/denominator
  print(paste("The value of ttc (ng/ml) is", round(result,2), sep=" "))
  return(round(result,2))
}

# dose rate
dr.fn <- function(ttc, cl, wtkg, f) {
  
  numerator   <- ttc * cl * wtkg
  denominator <- f
  result      <- (numerator/denominator) * 0.001
  print(paste("The value of dr (mg/h) is", round(result,3), sep=" "))
  return(round(result,3))
  
}

# dosing interval
di.fn <- function (msc, mec, ke) {
  numerator   <- log(msc/mec)
  denominator <- ke
  result      <- numerator/denominator
  print(paste("The value of di (h) is", round(result,1), sep=" "))
  return(round(result,1))
  
}

# back calculation TTC
bc.ttc.fn <- function (dr, f, cl, wtkg) {
  numerator   <- dr * 1000 * f
  denominator <- cl * wtkg
  result      <- numerator/denominator
  print(paste("The value of bc.ttc (mg/L) is", round(result*0.001,3), sep=" "))
  print(paste("The value of bc.ttc (ng/ml) is", round(result,3), sep=" "))
  return(round(result,2))
}

# accumulation rate
ar.fn <- function(ke, di) {
  numer  <- 1
  denom  <- 1 - exp(-ke * di)
  result <- numer/denom
  print(paste("The value of ar is", round(result,3), sep=" "))
  return(round(result,3))
}

# Oral dose
dpo.fn <- function(dr, di) {
  result      <- dr * di
  print(paste("The value of dpo (mg) is", round(result,2), sep=" "))
  print(paste("The value of dpo (ug) is", round(result*1000,3), sep=" "))
  return(round(result, 3))
}

# maintenance dose
dm.fn <- function (dr, di) {
  result <- dr * di
  print(paste("The value of dm (mg every di) is", round(result,2), sep=" "))
  print(paste("The value of dm (ug every di) is", round(result*1000,3), sep=" "))
  return(round(result,2))
}

# Cmax
cmax.fn <- function(f, dpo, vd, ar, wtkg) {
  numer <- f * dpo
  denom <- vd * wtkg
  result      <- (numer/denom) * ar
  print(paste("The value of cmax (mg/L) is", round(result/1000,3), sep=" "))
  print(paste("The value of cmax (ng/mL) is", round(result,2), sep=" "))
  return(round(result, 3))
}

# Cmin
cmin.fn <- function(cmax, ke, di) {
  result      <- cmax * exp(-ke * di) 
  print(paste("The value of cmin (mg/L) is", round(result/1000,3), sep=" "))
  print(paste("The value of cmin (ng/mL) is", round(result,2), sep=" "))
  return(round(result, 3))
}

# Css
css.fn <- function(f, dpo, di, cl, ar, wtkg) {
  numer  <- f * (dpo/di)
  denom  <- cl * wtkg
  result <- (numer/denom) * ar 
  print(paste("The value of css (mg/L) is", round(result/1000,3), sep=" "))
  print(paste("The value of css (ng/mL) is", round(result,2), sep=" "))
  return(round(result, 3))
}


# Ct - clearance curve (concentration time curve)
ct.fn <- function (cmax, ke, time=0) {
        result <- cmax * exp(-ke * time)
        print(paste("The value of ct (mg/L) at t h is", round(result/1000, 3), sep = " "))
		print(paste("The value of ct (ng/mL) at t h is", round(result, 2), sep = " "))
		return(round(result, 3))
}

# loading dose cmax
dlcmax.fn <- function(cmax,vd){
              result <- cmax * vd
              print(paste("The value of dlcmax (mg) is", round(result, 2), 
              sep = " "))
			  return(round(result, 2))

}

# loading dose ar
dlar.fn <- function(dm, ar){
            result <- dm * ar
            print(paste("The value of dlar (mg) is", round(result, 2), 
              sep = " "))
			return(round(result, 2))

}

# infusion rate R0 (mg/h)
R0.fn <- function(css, cl) {
           result <- css * cl
           print(paste("The value of R0 (ug/hr) is", round(result, 3), sep = " "))
		   return(round(result, 3))

}
####################### END ####################### 
