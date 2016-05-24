# File read.velocities.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigairini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                         'read_velocities'                                    # 
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates:                                                                     #        
################################################################################
 
                      
read_velocities <- function(file="Velocities.txt", ... ) {
                         
   out <- read_particles(file=file,...)   
   names(out) <- c("velocities", "gofs", "best.velocity", "best.gof")   
   return(out)
  
}  # 'read_velocities' END
