#########################################################################
#      $Log: Skernel.S,v $
# Revision 1.2  1995/04/05  18:56:55  bruno
# *** empty log message ***
#
# Revision 1.1  1995/04/02  01:04:16  bruno
# Initial revision
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################





rwkernel <- function(node, phinode, nvoice, x.inc = 1, x.min = node[1],
	x.max = node[length(node)], w0 = 2*pi, plot = FALSE)
#########################################################################
#     rwkernel:   
#     ------
#      Computes the cost from the sample of points on the estimated ridge
#      and the matrix used in the reconstruction of the original signal
#      The output is a lng x lng matrix of complex numbers, lng being the
#      number of points at which the signal is to be reconstructed. 
#      The dependence upon the original signal is only through the 
#      sample (node,phinode) of the ridge. This should be the output
#      of a previously run estimation procedure.
#
#     Input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the scale variable a for the nodes of the ridge
#      nvoice: number of scales within 1 octave
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      w0: central frequency of the wavelet
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     Output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{

cat("x.min = ",x.min,"\n")
cat("x.max = ",x.max,"\n")
cat("x.inc = ",x.inc,"\n")
   lng <- as.integer((x.max-x.min)/x.inc)+1
   nbnode <- length(node)

   cat("lng = ",lng,"\n");

   b.start <- x.min - 50
   b.end <- x.max + 50
   
   ker.r <- matrix(0, lng, lng)
   ker.i <- matrix(0, lng, lng)

   dim(ker.r) <- c(lng * lng, 1)
   dim(ker.i) <- c(lng * lng, 1)
   phinode <-  2 * 2^(phinode/nvoice)
  

  z <- .C("kernel",
        ker.r = as.double(ker.r),
        ker.i = as.double(ker.i),
        as.integer(x.min),
        as.integer(x.max),
	as.integer(x.inc),
	as.integer(lng),
        as.double(node),
        as.double(phinode),
        as.integer(nbnode),
        as.double(w0),
        as.double(b.start),
        as.double(b.end),
           PACKAGE="Rwave")


  ker.r <- z$ker.r
  ker.i <- z$ker.i
  dim(ker.r) <- c(lng, lng)
  dim(ker.i) <- c(lng, lng)
  ker <- matrix(0, lng, lng)
  ker <- ker.r + 1i*ker.i

  if(plot == TRUE){
     par(mfrow=c(1,1))
     image(Mod(ker))
     title("Matrix of the reconstructing kernel (modulus)")
  }

  ker
}

rkernel <- function(node, phinode, nvoice, x.inc = 1, x.min = node[1],
	x.max = node[length(node)], w0 = 2*pi, plot = FALSE)
#########################################################################
#     rkernel:   
#     -------
#      Same as kernel, except that a real valued kernel is computed
#      (precisely the real part of the previous kernel); this applies to
#      real signals.
#
#     Input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the scale variable a for the nodes of the ridge
#      nvoice: number of scales within 1 octave
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      w0: central frequency of the wavelet
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     Output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{
   lng <- as.integer((x.max-x.min)/x.inc)+1
   nbnode <- length(node)

   b.start <- x.min - 50
   b.end <- x.max + 50

   ker <- matrix(0, lng, lng)
   dim(ker) <- c(lng * lng, 1)
   phinode <-  2 * 2^(phinode/nvoice)

  z <- .C("rkernel",
        ker = as.double(ker),
        as.integer(x.min),
        as.integer(x.max),
	as.integer(x.inc),
	as.integer(lng),
        as.double(node),
        as.double(phinode),
        as.integer(nbnode),
        as.double(w0),
        as.double(b.start),
        as.double(b.end),
           PACKAGE="Rwave")

  ker <- z$ker
  dim(ker) <- c(lng, lng)

  if(plot == TRUE){
     par(mfrow=c(1,1))
     image(Mod(ker))
     title("Matrix of the Q2 kernel (modulus)")
  }

  ker
}



fastkernel <- function(node, phinode, nvoice, x.inc = 1, x.min = node[1],
	x.max = node[length(node)], w0 = 2*pi, plot = FALSE)
#########################################################################
#     fastkernel:   
#     -----------
#	   Same as kernel, except that the kernel is computed
#	     using Riemann sums instead of Romberg integration.
#
#     Input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the scale variable a for the nodes of the ridge
#      nvoice: number of scales within 1 octave
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      w0: central frequency of the wavelet
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     Output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{
   lng <- as.integer((x.max-x.min)/x.inc)+1
   nbnode <- length(node)

   b.start <- x.min - 50
   b.end <- x.max + 50

   ker.r <- matrix(0, lng, lng)
   ker.i <- matrix(0, lng, lng)
   dim(ker.r) <- c(lng * lng, 1)
   dim(ker.i) <- c(lng * lng, 1)
   phinode <-  2 * 2^(phinode/nvoice)

  z <- .C("fastkernel",
        ker.r = as.double(ker.r),
        ker.i = as.double(ker.i),
        as.integer(x.min),
        as.integer(x.max),
	as.integer(x.inc),
	as.integer(lng),
        as.double(node),
        as.double(phinode),
        as.integer(nbnode),
        as.double(w0),
        as.double(b.start),
        as.double(b.end),
           PACKAGE="Rwave")


  ker.r <- z$ker.r
  ker.i <- z$ker.i
  dim(ker.r) <- c(lng, lng)
  dim(ker.i) <- c(lng, lng)
  ker <- matrix(0, lng, lng)
  ker <- ker.r + 1i*ker.i

  if(plot == TRUE){
     par(mfrow=c(1,1))
     image(Mod(ker))
     title("Matrix of the reconstructing kernel (modulus)")
  }

  ker
}



zerokernel <- function(x.inc = 1, x.min,x.max)
#########################################################################
#     zerokernel:   
#     -----------
#	Returns a zero kernel.
#
#     Input:
#     ------
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#
#     Output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{
   lng <- as.integer((x.max-x.min)/x.inc)+1

   ker <- matrix(0, lng, lng)

   ker
}




gkernel <- function(node, phinode, freqstep, scale, x.inc = 1,
	x.min = node[1], x.max = node[length(node)], plot = FALSE)
#########################################################################
#     gkernel:   
#     -------
#	Same as kernel, for the case of Gabor transform.
#
#     input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the frequency variable for the ridge nodes
#      freqstep: sampling rate for the frequency axis
#      scale: size of the window
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{
   lng <- as.integer((x.max-x.min)/x.inc)+1
   nbnode <- length(node)
   b.start <- node[1] - 50
   b.end <- node[length(node)] + 50
   ker <- matrix(0, lng, lng)
   dim(ker) <- c(lng * lng, 1)

   phinode <- phinode * freqstep


  z <- .C("gkernel",
        ker = as.double(ker),
        as.integer(x.min),
        as.integer(x.max),
	as.integer(x.inc),
	as.integer(lng),
        as.double(node),
        as.double(phinode),
        as.integer(nbnode),
        as.double(scale),
        as.double(b.start),
        as.double(b.end),
           PACKAGE="Rwave")


  ker <- z$ker
  dim(ker) <- c(lng, lng)

  if(plot == TRUE){
     par(mfrow=c(1,1))
     image(Mod(ker))
     title("Matrix of the reconstructing kernel (modulus)")
  }

  ker
}



fastgkernel <- function(node, phinode, freqstep, scale, x.inc = 1,
	x.min = node[1], x.max = node[length(node)], plot = FALSE)
#########################################################################
#     fastgkernel:   
#     ------------
#	Same as gkernel, except that the kernel is computed
#	  using Riemann sums instead of Romberg integration.
#
#
#     input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the frequency variable for the ridge nodes
#      freqstep: sampling rate for the frequency axis
#      scale: size of the window
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
{
   lng <- as.integer((x.max-x.min)/x.inc)+1
   nbnode <- length(node)

   b.start <- x.min - 50
   b.end <- x.max + 50

   ker <- matrix(0, lng, lng)
   dim(ker) <- c(lng * lng, 1)
   phinode <-  phinode*freqstep


  z <- .C("fastgkernel",
        ker = as.double(ker),
        as.integer(x.min),
        as.integer(x.max),
	as.integer(x.inc),
	as.integer(lng),
        as.double(node),
        as.double(phinode),
        as.integer(nbnode),
        as.double(scale),
        as.double(b.start),
        as.double(b.end),
           PACKAGE="Rwave")


  ker <- z$ker
  dim(ker) <- c(lng, lng)

  if(plot == TRUE){
     par(mfrow=c(1,1))
     image(Mod(ker))
     title("Matrix of the reconstructing kernel (modulus)")
  }

  ker
}









