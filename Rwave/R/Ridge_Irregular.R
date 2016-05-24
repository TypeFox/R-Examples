#########################################################################
#      $Log: Ridge_Recons.S,v $
# Revision 1.2  1995/04/05  18:56:55  bruno
# *** empty log message ***
#
# Revision 1.1  1995/04/02  01:04:16  bruno
# Initial revision
#
#
#               (c) Copyright  1997                             
#                          by     
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





irregrec <- function(siginput, cwtinput,phi,compr,noct,nvoice,
                     epsilon = 0.5, w0 = 2*pi, prob = 0.8, fast = FALSE,
                     plot = FALSE, para = 0, hflag = FALSE, check = FALSE,
                     minnbnodes = 2, real = FALSE)
#########################################################################
#     irregrec:
#     ---------
#      Reconstruction of a real valued signal from a (continuous) ridge
#      (uses a irregular sampling of the ridge)
#
#      input:
#      -------
#      siginput: input signal
#      cwtinput: Continuous wavelet transform (output of cwt)
#      phi: (unsampled) ridge
#      compr: subsampling rate for the wavelet coefficients (at scale 1)
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      epsilon: coeff of the Q2 term in reconstruction kernel
#      w0: central frequency of Morlet wavelet
#      fast: if set to TRUE, the kernel is computed using Riemann
#           sums instead of Romberg's quadrature
#      plot: if set to TRUE, displays original and reconstructed signals
#      para: constant for extending the support of reconstructed signal
#            outside the ridge. 
#      check: if set to TRUE, computes the wavelet transform of the
#            reconstructed signal
#      minnbnodes: minimum number of nodes for the reconstruction
#      real: if set to TRUE, only uses constraints on the real part
#            of the transform for the reconstruction.
#      prob: only keep the lambda's values greater than this percential after regular sampling
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#      solskel: wavelet transform of sol, restricted to the ridge
#      inputskel: wavelet transform of signal, restricted to the ridge
#      Q2: second part of the reconstruction kernel
#      nbnodes: number of nodes used for the reconstruction.
#
##########################################################################
{

#
#  Generate (regularly) sampled ridge
#


   tmp <- wRidgeSampling(phi,compr,nvoice)

   node <- tmp$node
cat("node at = ", node, "\n")
   phinode <- tmp$phinode
   nbnodes <- tmp$nbnodes
   phinode <- as.integer(phinode)	

   if(nbnodes < minnbnodes){
      cat(" Chain too small\n")
      NULL 
   }

   phi.x.min <- 2 * 2^(phinode[1]/nvoice)
   phi.x.max <- 2 * 2^(phinode[length(node)]/nvoice)

   x.min <- node[1]
   x.max <- node[length(node)]
   x.max <- x.max + round(para * phi.x.max)
   x.min <- x.min - round(para * phi.x.min)

   node <- node + 1 - x.min
   x.inc <- 1
   np <- as.integer((x.max - x.min)/x.inc) +1

    # Generating the Q2 term in reconstruction kernel
    if(epsilon == 0)
       Q2 <- 0
    else {
       if (fast == FALSE)
         Q2 <- rkernel(node, phinode, nvoice, x.min = x.min,
                       x.max = x.max,w0 = w0)
       else
	 Q2 <- fastkernel(node, phinode, nvoice, x.min = x.min,
                          x.max = x.max,w0 = w0)
     }
    cat(" kernel; ")

    # Generating the Q1 term in reconstruction kernel
    if (hflag == TRUE){
       one <- numeric(np)
       one[] <- 1
     }
    else{
       one <- numeric(np)
       one[] <- 1
    }

    if (epsilon !=0 ){
      Q <-  epsilon * Q2
      for(j in 1:np)
         Q[j,j] <- Q[j,j] + one[j]
      Qinv <- solve(Q)
    }
    else{
       Qinv <-  1/one
     }

    tmp2 <- ridrec(cwtinput, node, phinode, noct, nvoice, Qinv,
                   epsilon, np, w0 = w0, check = check, real = real)

# 
# Now perform the irregular sampling 
#
# C and Splus ...

   tmp <- RidgeIrregSampling(phi, node, tmp2$lam, prob)
   node <- tmp$node

  cat("node at  ", node,"\n")

   phinode <- tmp$phinode

  cat("phinode at  ", phinode,"\n")
   nbnodes <- tmp$nbnodes
   phinode <- as.integer(phinode)	

   if(nbnodes < minnbnodes){
      cat(" Chain too small\n")
      NULL 
   }

   phi.x.min <- 2 * 2^(phinode[1]/nvoice)
   phi.x.max <- 2 * 2^(phinode[length(node)]/nvoice)

   x.min <- node[1]
   x.max <- node[length(node)]

   x.max <- x.max + round(para * phi.x.max)
   x.min <- x.min - round(para * phi.x.min)

   node <- node + 1 - x.min
   x.inc <- 1


   np <- as.integer((x.max - x.min)/x.inc) +1
   cat("(size:",np,",",nbnodes,"nodes):")

    # Generating the Q2 term in reconstruction kernel
    if(epsilon == 0)
       Q2 <- 0
    else {
       if (fast == FALSE)
         Q2 <- rkernel(node, phinode, nvoice, x.min = x.min,
                       x.max = x.max,w0 = w0)
       else
	 Q2 <- fastkernel(node, phinode, nvoice, x.min = x.min,
                          x.max = x.max,w0 = w0)
     }
    cat(" kernel; ")

    # Generating the Q1 term in reconstruction kernel
    if (hflag == TRUE){
       one <- numeric(np)
       one[] <- 1
     }
    else{
       one <- numeric(np)
       one[] <- 1
    }

    if (epsilon !=0 ){
      Q <-  epsilon * Q2
      for(j in 1:np)
        Q[j,j] <- Q[j,j] + one[j]
      Qinv <- solve(Q)
    }
    else{
      Qinv <-  1/one
    }


    tmp2 <- ridrec(cwtinput, node, phinode, noct, nvoice, Qinv,
                   epsilon, np, w0 = w0, check = check, real = real)


    if(plot == TRUE){
       par(mfrow=c(2,1))
       plot.ts(Re(siginput))
       title("Original signal")
       plot.ts(Re(tmp2$sol))
       title("Reconstructed signal")
    }


    list(sol = tmp2$sol, A = tmp2$A, lam = tmp2$lam, dualwave = tmp2$dualwave,
	  morvelets = tmp2$morvelets, solskel = tmp2$solskel,
          inputskel = tmp2$inputskel, Q2 = Q2, nbnodes = nbnodes)
}



