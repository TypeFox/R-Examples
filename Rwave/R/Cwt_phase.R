#########################################################################
#      $Log: Cwt_phase.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################





cwtp <- function(input, noctave, nvoice = 1, w0 = 2*pi,
	twoD = TRUE, plot = TRUE)
#########################################################################
#       cwtp:   
#       ----
#        continuous wavelet transform function:
# 	  compute the continuous wavelet transform with Morlet wavelet,
#    	  and its phase derivative (minus the wavelet frequency)
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 noctave: number of powers of 2 for the scale variable
#	 nvoice: number of scales between 2 consecutive powers of 2
#        w0: central frequency of Morlet wavelet
#	 twoD: if set to TRUE, organizes the output as a 2D array 
#			(signal_size X nb_scales)
#		      if not: 3D array (signal_size X noctave X nvoice)
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       output:
#       -------
#        wt: continuous (complex) wavelet transform
#        f: derivative of wavelet transform phase
#
#########################################################################
{
   oldinput <- input
   isize <- length(oldinput)

   tmp <- adjust.length(oldinput)
   input <- tmp$signal
   newsize <- length(input)

   pp <- noctave * nvoice
   Routput <- matrix(0,newsize,pp)
   Ioutput <- matrix(0,newsize,pp)
   Foutput <- matrix(0,newsize,pp)
   output <- matrix(0,newsize,pp)
   dim(Routput) <- c(pp * newsize,1)
   dim(Ioutput) <- c(pp * newsize,1)
   dim(Foutput) <- c(pp * newsize,1)
   dim(input) <- c(newsize,1)

   ## dyn.load("/usr/people/bruno/Swave/wave.a")
   ## dyn.load("/usr/people/whwang/Swave/wave.a")

   z <- .C("Scwt_phase",
           as.double(input),
           Rtmp = as.double(Routput),
           Itmp = as.double(Ioutput),
           Ftmp = as.double(Foutput),
           as.integer(noctave),
           as.integer(nvoice),
           as.integer(newsize),
           as.double(w0),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   Foutput <- z$Ftmp
   dim(Routput) <- c(newsize,pp)
   dim(Foutput) <- c(newsize,pp)
   dim(Ioutput) <- c(newsize,pp)

   if(twoD) {
     Ftmp <- Foutput[1:isize,]
     output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
     if(plot) image(Mod(output))
     list(wt=output,f=Ftmp)
   } 
   else {
     Rtmp <- array(0,c(isize,noctave,nvoice))
     Itmp <- array(0,c(isize,noctave,nvoice))
     Ftmp <- array(0,c(isize,noctave,nvoice))
     for(i in 1:noctave)
       for(j in 1:nvoice) {
         Rtmp[,i,j] <- Routput[1:isize,(i-1)*nvoice+j]
         Itmp[,i,j] <- Ioutput[1:isize,(i-1)*nvoice+j]
         Ftmp[,i,j] <- Foutput[1:isize,(i-1)*nvoice+j]
         }       
     list(wt=Rtmp+1i*Itmp, f=Ftmp)
   }
}





