#########################################################################
#      $Log: Cwt_Squeezing.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################






cwtsquiz <- function(input, noctave, nvoice = 1, w0 = 2*pi,
	twoD = TRUE, plot = TRUE)
#########################################################################
#      cwtsquiz:
#      --------
#       squeezed wavelet transform function
#
#       Input:
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
#        tmp: continuous (complex) squeezed wavelet transform
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
   output <- matrix(0,newsize,pp)
   dim(Routput) <- c(pp * newsize,1)
   dim(Ioutput) <- c(pp * newsize,1)
   dim(input) <- c(newsize,1)

   z <- .C("Scwt_squeezed",
           as.double(input),
           Rtmp = as.double(Routput),
           Itmp = as.double(Ioutput),
           as.integer(noctave),
           as.integer(nvoice),
           as.integer(newsize),
           as.double(w0),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   dim(Routput) <- c(newsize,pp)
   dim(Ioutput) <- c(newsize,pp)
   if(twoD) {
     output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
     if(plot){
        image(Mod(output),xlab="Time", ylab="log(scale)")
	title("Squeezed Wavelet Transform Modulus")
     }
     output
   } 
   else {
     Rtmp <- array(0,c(isize,noctave,nvoice))
     Itmp <- array(0,c(isize,noctave,nvoice))
     for(i in 1:noctave)
       for(j in 1:nvoice) {
         Rtmp[,i,j] <- Routput[1:isize,(i-1)*nvoice+j]
         Itmp[,i,j] <- Ioutput[1:isize,(i-1)*nvoice+j]
      }
     Rtmp + 1i*Itmp
   }
}






