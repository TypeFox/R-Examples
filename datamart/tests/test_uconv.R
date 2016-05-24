library(datamart)

test_uconv <- function() {
  
  # simple use cases
  if(abs(uconv(1, "TWh", "PJ", "Energy")-3.6)>0.0001) cat("Energy Unitset failed.\n")
  if(abs(uconv(13, "km", "m", "Length")-13000)>0.0001) cat("Length Unitset failed.\n")
  if(abs(uconv(-5, "t", "kg", "Mass")+5000)>0.0001) cat("Mass Unitset failed.\n")
  if(abs(uconv(100*100, "m\u00b2", "ha", "Area")-1)>0.0001) cat("Area Unitset failed.\n")
  if(abs(uconv(1, "l", "dm\u00b3", "Volume")-1)>0.0001) cat("Volume Unitset failed.\n")
  
  # vector use cases
  if(any(abs(uconv(1:3, "TWh", "PJ", "Energy")-((1:3)*3.6))>0.0001)) cat("Vector Energy Unitset failed.\n")
  if(any(abs(uconv(c(1313, 13), "km", "m", "Length")-(c(1313, 13)*1000))>0.0001)) cat("Vector Length Unitset failed.\n")
  if(any(abs(uconv(c(-5,4), "t", "kg", "Mass")-(c(-5,4)*1000))>0.0001)) cat("Vector Mass Unitset failed.\n")
  if(any(abs(uconv(100*100*c(1,3,5), "m\u00b2", "ha", "Area")-c(1,3,5))>0.0001)) cat("Vector Area Unitset failed.\n")
  if(any(abs(uconv(1:3, "l", "dm\u00b3", "Volume")-(1:3))>0.0001)) cat("Vector Volume Unitset failed.\n")
  
  # guess unitset
  w_opt <- getOption("warn")
  options(warn=-1)
  if(abs(uconv(1, "TWh", "PJ")-3.6)>0.0001) cat("Guess Energy Unitset failed.\n")
  if(abs(uconv(13, "km", "m")-13000)>0.0001) cat("Guess Length Unitset failed.\n")
  if(abs(uconv(-5, "t", "kg")+5000)>0.0001) cat("Guess Mass Unitset failed.\n")
  if(abs(uconv(100*100, "m\u00b2", "ha")-1)>0.0001) cat("Guess Area Unitset failed.\n")
  if(abs(uconv(1, "l", "dm\u00b3")-1)>0.0001) cat("Guess Volume Unitset failed.\n")
  options(warn=w_opt)
  
  # NAs
  res <- uconv(c(1,NA,3), "TWh", "PJ", "Energy")
  if(any(abs(res[c(1,3)]-(c(1,3)*3.6))>0.0001) || !is.na(res[2])) cat("NA uconv failed.\n")
  
  
  # finish
  cat("Done.\n")
}

test_uconv()
