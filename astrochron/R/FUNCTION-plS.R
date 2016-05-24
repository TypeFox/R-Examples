### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### plS: set up plotting parameters for vertical stratigraphic plots. This 
###      is ususally invoked after function pl (SRM: January 15, 2014)
###
###########################################################################

plS <- function (f=T,s=1)
{

if(f) cat("\n* Setting parameters for vertical stratigraphic data plots- First plot.\n")
if(!f) cat("\n* Setting parameters for vertical stratigraphic data plots- Remaining plots.\n")

if(f) par(mar = c(5,5,5,0.2))

# set font and symbol size
   par(cex=s)
   par(cex.axis=s)
   par(cex.lab=s)
   par(cex.main=s)

# do not plot box
par(bty="n")

if(!f) 
 {
   par(mar = c(5,0.2,5,0.2))
# do not annotate y-axis
   par(yaxt="n")
 }  

### END function plS
}
