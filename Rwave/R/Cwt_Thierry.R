#########################################################################
#       $Log: Cwt_Thierry.S,v $
#########################################################################
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################





cwtTh <- function(input, noctave, nvoice = 1, moments,
	twoD = TRUE, plot = TRUE)
#########################################################################
#       cwtTh:   Cauchy's wavelet transform
#       -----
# 	 continuous wavelet transform function:
#         compute the continuous wavelet transform with (complex-valued)
#	  Cauchy's wavelet
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 noctave: number of powers of 2 for the scale variable
#	 nvoice: number of scales between 2 consecutive powers of 2
#        moments: number of vanishing moments for the wavelet
#	 twoD: if set to TRUE, organizes the output as a 2D array 
#			(signal_size X nb_scales)
#		      if not: 3D array (signal_size X noctave X nvoice)
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       output:
#       -------
#        tmp: continuous (complex) wavelet transform
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

   z <- .C("Scwt_thierry",
           as.double(Re(input)),
           as.double(Im(input)),
           Rtmp = as.double(Routput),
           Itmp = as.double(Ioutput),
           as.integer(noctave),
           as.integer(nvoice),
           as.integer(newsize),
           as.integer(moments),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   dim(Routput) <- c(newsize,pp)
   dim(Ioutput) <- c(newsize,pp)
   if(twoD) {
     output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
     if(plot) { 
       image(Mod(output),xlab="Time",ylab="log(scale)")
       title("Wavelet Transform Modulus")
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



vwtTh <- function(input, scale, moments)
#########################################################################
#       vwtTh:   
#       -----
#        continuous wavelet transform on one scale:
# 	 compute the continuous wavelet transform with (complex-valued)
#	  Cauchy's wavelet
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#        scale: value of the scale at which the transform is computed
#        moments: number of vanishing moments
#
#       output:
#       -------
#        Routput + i Ioutput: voice wavelet transform (complex 1D array)
#
#########################################################################
{
   oldinput <- input
   isize <- length(oldinput)

   tmp <- adjust.length(oldinput)
   input <- tmp$signal
   newsize <- length(input)

   Routput <- numeric(newsize)
   Ioutput <- numeric(newsize)
   dim(input) <- c(newsize,1)

   z <- .C("Svwt_Thierry",
           as.double(Re(input)),
           as.double(Im(input)),
           Rtmp = as.double(Routput),
           Itmp = as.double(Ioutput),
           as.double(scale),
           as.integer(newsize),
           as.integer(moments),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   Routput[1:isize] + 1i*Ioutput[1:isize]
}



