### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### getLaskar: download Laskar et al. (2004; 2011) astronomical solutions. 
###          (SRM: January 22, 2015)
###
###########################################################################


getLaskar <- function (sol="la04")
{
        i=0
        if(sol=="la04" || sol=="la10a" || sol=="la10b" || sol=="la10c" || sol=="la10d" || sol=="la11") i=1
        if(i==0) 
         {
            cat("\n**** ERROR: specified solution is not available.\n")
            stop("    TERMINATING NOW!")
         }
        
        tempLaskar <- tempfile()
# download data
        if(sol=="la04")
        {
          cat(" * Downloading Laskar et al. (2004) astronomical solution\n\n")
          cat("   Please cite: Laskar, J., Robutel, P., Joutel, F., Gastineau, M.,\n") 
          cat("   Correia, A.C.M., Levrard, B., 2004, A long term numerical solution \n")
          cat("   for the insolation quantities of the Earth: Astron. Astrophys., Volume 428, 261-285.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la04.txt.bz2",tempLaskar)
        }
  
        if(sol=="la10a")
        {
          cat(" * Downloading Laskar et al. (2011) astronomical solution: 10a\n\n")
          cat("   Please cite: Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011,\n") 
          cat("   La2010: A new orbital solution for the long-term motion of the Earth:\n")
          cat("   Astron. Astrophys., Volume 532, A89.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la10a.txt.bz2",tempLaskar)
        }

        if(sol=="la10b")
        {
          cat(" * Downloading Laskar et al. (2011) astronomical solution: 10b\n\n")
          cat("   Please cite: Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011,\n") 
          cat("   La2010: A new orbital solution for the long-term motion of the Earth:\n")
          cat("   Astron. Astrophys., Volume 532, A89.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la10b.txt.bz2",tempLaskar)
        }

        if(sol=="la10c")
        {
          cat(" * Downloading Laskar et al. (2011) astronomical solution: 10c\n\n")
          cat("   Please cite: Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011,\n") 
          cat("   La2010: A new orbital solution for the long-term motion of the Earth:\n")
          cat("   Astron. Astrophys., Volume 532, A89.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la10c.txt.bz2",tempLaskar)
        }

        if(sol=="la10d")
        {
          cat(" * Downloading Laskar et al. (2011) astronomical solution: 10d\n\n")
          cat("   Please cite: Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011,\n") 
          cat("   La2010: A new orbital solution for the long-term motion of the Earth:\n")
          cat("   Astron. Astrophys., Volume 532, A89.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la10d.txt.bz2",tempLaskar)
        }

        if(sol=="la11")
        {
          cat(" * Downloading Laskar et al. (2011) astronomical solution: 11\n\n")
          cat("   Please cite: Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011,\n") 
          cat("   La2010: A new orbital solution for the long-term motion of the Earth:\n")
          cat("   Astron. Astrophys., Volume 532, A89.\n")
          cat("  AND:\n")
          cat("   Laskar, J., Gastineau, M., Delisle, J.-B., Farres, A., Fienga, A.: 2011,\n")
          cat("   Strong chaos induced by close encounters with Ceres and Vesta:\n")
          cat("   Astron. Astrophys., Volume 532, L4.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la11.txt.bz2",tempLaskar)
        }

        cat(" * Decompressing solution\n")
        la <- read.table(bzfile(tempLaskar),header=T)
        unlink(tempLaskar)
        if(sol != "la11") la=data.frame(cbind((1:249001)-1,la))
        if(sol == "la11") la=data.frame(cbind((1:100000)-1,la))
        colnames(la)[1] = c("Time_ka")
  
        return(la)

### END function getLaskar
}
