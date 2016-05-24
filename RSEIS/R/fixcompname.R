`fixcompname` <-
function(comp)
  {
##X## convert wierd component names to something more useful for seismic

    tcomp = "XXX"
    if(comp=="SHV"|| comp=="4"|| comp=="1" || comp=="V" || comp=="v" || comp=="G1V") { return("Vertical") }
    if(comp=="SHN"|| comp=="5"|| comp=="2" || comp=="N" || comp=="n" || comp=="G1N") { return("North") }
    if(comp=="SHE"|| comp=="6"|| comp=="3" || comp=="E" || comp=="e" || comp=="G1E") { return("East") }

    

 return(tcomp)
  }

