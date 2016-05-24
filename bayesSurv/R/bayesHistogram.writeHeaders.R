#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       bayesHistogram.writeHeaders.R       ####
####                                                 ####
#### FUNCTIONS:  bayesHistogram.writeHeaders         ####
#########################################################

### ======================================
### bayesHistogram.writeHeaders
### ======================================
## Subfunction for bayesHistogram.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayesHistogram.writeHeaders <- function(dir, design, prior.init, store)
{
   nP <- length(design$status)/design$dim
  
   sink(paste(dir, "/iteration.sim", sep = ""), append = FALSE)
   cat("iteration", "\n"); sink()

   write.headers.Gspline(dir=dir, dim=design$dim, nP=nP,
                         label="", gparmi=prior.init$Gparmi, store.a=store$a, store.y=store$y, store.r=store$r, care.of.y=TRUE)
}  
