### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function getcolor - (SRM: August 2, 2012, August 6, 2013)
###
### query R for color information
###########################################################################

getColor <- function (color) 
{

# names
  out1 = colors()[grep(color,colors())]
# IDs
  out2 = grep(color,colors())

  cat("\nCOLORS CONTAINING THE WORD",color,":\n")
  out3=data.frame(cbind(out1,out2))
  colnames(out3)<- c('name','id')
  print(out3)
  
### END function getColor
}
