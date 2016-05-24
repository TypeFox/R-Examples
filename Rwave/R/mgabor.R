#########################################################################
#      $Log: MGabor.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





mcgt <- function(input, nvoice, freqstep = (1/nvoice), nscales = 10,
	scalestep = 5, initscale = 0, crit = 0, plot = TRUE,
	tchatche = FALSE)
#########################################################################
#       mcgt:   
#       ----
#        multiwindow continuous Gabor transform function:
# 	  compute the continuous Gabor transform with gaussian windows.
#	  selects the optimal window using L^p norm or entropy
#         optimization
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 nvoice: number of frequency bands
#        freqstep: sampling rate for the frequency axis
#        nscales: number of scales considered
#        scalestep: increment for scale parameter
#        initscale: initial value form the scale parameter
#        crit: criterion for optimization (L^p norm or entropy)
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#        rchatche: if set to TRUE, prints intermediate results
#
#       output:
#       -------
#        output: optimal continuous (complex) Gabor transform
#
#########################################################################
{
   oldinput <- input
   isize <- length(oldinput)

   tmp <- adjust.length(oldinput)
   input <- tmp$signal
   newsize <- length(input)

   pp <- nvoice
   Routput <- matrix(0,newsize,pp)
   Ioutput <- matrix(0,newsize,pp)
   output <- matrix(0,newsize,pp)
   dim(Routput) <- c(pp * newsize,1)
   dim(Ioutput) <- c(pp * newsize,1)
   dim(input) <- c(newsize,1)


   norm <- 1
   lpoptnorm <- 0
   entoptnorm <- 100000000

   optsca <- 0


#####################
# Loop over scales: #
#####################

   sca <- initscale
   for(k in 1:nscales){

      sca <- sca + scalestep

# Compute Gabor transform
# -----------------------
      z <- .C("Sgabor",
              as.double(input),
              Rtmp = as.double(Routput),
              Itmp = as.double(Ioutput),
              as.integer(nvoice),
              as.double(freqstep),
              as.integer(newsize),
              as.double(sca),
              PACKAGE="Rwave")

      Routput <- z$Rtmp
      Ioutput <- z$Itmp
      dim(Routput) <- c(newsize,pp)
      dim(Ioutput) <- c(newsize,pp)

# Compute L^2 norm for normalization
# ----------------------------------
      pexp <- as.double(2)
      z <- .C("Lpnorm",
              ltwonorm = as.double(norm),
              as.double(pexp),
              as.double(Routput),
              as.double(Ioutput),
              as.integer(newsize),
              as.integer(nvoice),
              PACKAGE="Rwave")
      ## cat("l2 norm=",z$ltwonorm,"\n")

# Normalize
# ---------
      Routput <- Routput/z$ltwonorm
      Ioutput <- Ioutput/z$ltwonorm

      if(crit == 0) {
         z <- .C("entropy",
                 lpnorm = as.double(norm),
                 as.double(Routput),
                 as.double(Ioutput),
                 as.integer(newsize),
                 as.integer(nvoice),
                 PACKAGE="Rwave")
         if(tchatche) {
            cat("     scale=", sca, "; entropy=", z$lpnorm, "\n")
         }

         if (z$lpnorm < entoptnorm){
            entoptnorm <- z$lpnorm
            optsca <- sca
            output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
         }
      optnorm <- entoptnorm
      }

      else {

        ## Compute L^p norm
        ## ----------------
        pexp <- as.double(crit)
        z <- .C("Lpnorm",
                lpnorm = as.double(norm),
                as.double(pexp),
                as.double(Routput),
                as.double(Ioutput),
                as.integer(newsize),
                as.integer(nvoice),
                PACKAGE="Rwave")
        if(tchatche) {
          cat("     scale=", sca,"; l", pexp," norm=", z$lpnorm, "\n")
        }
        
        if(z$lpnorm > lpoptnorm) {
          lpoptnorm <- z$lpnorm
          optsca <- sca
          output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
        }
        optnorm <- lpoptnorm
      }
    }
   
   cat("   Optimal scale: ", optsca, "\n")
   
   if(plot) {
     image(Mod(output), xlab="Time", ylab="Frequency")
     title("Gabor Transform Modulus")
   }
   
   output
}
