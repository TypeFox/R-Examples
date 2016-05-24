### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### getData: download data from astrochon server. 
###          (SRM: October 21, 2015)
###
###########################################################################


getData <- function (dat="1262-a*")
{
        
        tempDat <- tempfile()
# download data
        if(dat=="926B-18O")
        {
          cat(" * Downloading adjusted benthic foraminifera oxygen isotope data from ODP Site 926B\n\n")
          cat("   Please cite: Paelike, H., Frazier, J., Zachos, J.C., 2006,\n") 
          cat("   Extended orbitally forced palaeoclimatic records from the  \n")
          cat("   equatorial Atlantic Ceara Rise: Quaternary Science Reviews, 25(23-24), 3138-3149\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/926B-18O.txt.bz2",tempDat)
        }
  
        if(dat=="1262-a*")
        {
          cat(" * Downloading a* color data from ODP Site 1262\n\n")
          cat("   Please cite: Zachos, J.C., Kroon, D., Blum, P., et al., 2004,\n") 
          cat("   Proc. ODP, Init. Repts., 208: College Station, TX (Ocean Drilling\n")
          cat("   Program). doi:10.2973/odp.proc.ir.208.2004\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/1262_a.txt.bz2",tempDat)
        }

        cat(" * Decompressing\n")
        dat <- read.table(bzfile(tempDat),header=T)
  
        return(dat)

### END function getData
}
