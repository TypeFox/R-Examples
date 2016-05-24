#########################################################################
#      $Log: gRidge_Recons.S,v $
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
#                  University of California, Irvine             
#                  All right reserved                           
#########################################################################




girregrec <- function(siginput, gtinput, phi, nbnodes, nvoice,
	freqstep, scale, epsilon = 0.5,fast = FALSE, prob= 0.8, plot = FALSE,
	para = 0, hflag = FALSE, real = FALSE, check = FALSE)
#########################################################################
#     girregrec:
#     ---------
#      Reconstruction of a real valued signal from a (continuous) 
#      Gabor ridge (uses a irregular sampling of the ridge)
#
#      input:
#      -------
#      siginput: input signal
#      gtinput: Continuous gabor transform (output of cgt)
#      phi: (unsampled) ridge
#      nbnodes: number of nodes used for the reconstruction.
#      nvoice: number of different scales per octave
#      freqstep: sampling rate for the frequency axis
#      epsilon: coeff of the Q2 term in reconstruction kernel
#      fast: if set to TRUE, the kernel is computed using Riemann
#           sums instead of Romberg's quadrature
#      plot: if set to TRUE, displays original and reconstructed signals
#      para: no comment !
#      prob: sampling (irregularly) at the places of the ridge where their lambdas 
#            are more than prob in lambda distribution.
#      hflag: 
#      real: if set to TRUE, only uses constraints on the real part
#            of the transform for the reconstruction.
#      check: if set to TRUE, computes the wavelet transform of the
#            reconstructed signal
#      minnbnodes: minimum number of nodes for the reconstruction
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
#########################################################################
{

#
# first performing regular sampling
#

   tmp <- RidgeSampling(phi,nbnodes)

   node <- tmp$node
   phinode <- tmp$phinode

   phi.x.min <- scale
   phi.x.max <- scale
   
   x.min <- node[1]
   x.max <- node[length(node)]
   x.max <- x.max + round(para * phi.x.max)
   x.min <- x.min - round(para * phi.x.min)

   node <- node - x.min + 1
   
   x.inc <- 1
   np <- as.integer((x.max-x.min)/x.inc)+1

   cat(" (np:",np,")")

   if(epsilon == 0)
       Q2 <- 0
   else {
      if (fast == FALSE)
	Q2 <- gkernel(node,phinode,freqstep,scale,x.min = x.min, x.max = x.max)
      else
	Q2 <- fastgkernel(node,phinode,freqstep,scale,x.min = x.min,
		x.max = x.max)
      }
   cat(" kernel;")

# Generating the Q1 term in reconstruction kernel
   if (hflag == TRUE)
      one <- gsampleOne(node,scale,np)
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

   tmp2 <- gridrec(gtinput,node,phinode,nvoice,
	freqstep,scale,Qinv,epsilon,np, real = real, check = check)

# 
# Now perform the irregular sampling 
#  

   tmp <- RidgeIrregSampling(phinode, node, tmp2$lam, prob)

   node <- tmp$node
   phinode <- tmp$phinode

   phi.x.min <- scale
   phi.x.max <- scale
   
   x.min <- node[1]
   x.max <- node[length(node)]
   x.max <- x.max + round(para * phi.x.max)
   x.min <- x.min - round(para * phi.x.min)

   node <- node - x.min + 1
   
   x.inc <- 1
   np <- as.integer((x.max-x.min)/x.inc)+1

   cat(" (np:",np,")")

   if(epsilon == 0)
       Q2 <- 0
   else {
      if (fast == FALSE)
	Q2 <- gkernel(node,phinode,freqstep,scale,x.min = x.min, x.max = x.max)
      else
	Q2 <- fastgkernel(node,phinode,freqstep,scale,x.min = x.min,
		x.max = x.max)
      }
   cat(" kernel;")

# Generating the Q1 term in reconstruction kernel
   if (hflag == TRUE)
      one <- gsampleOne(node,scale,np)
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

   tmp2 <- gridrec(gtinput,node,phinode,nvoice,
	freqstep,scale,Qinv,epsilon,np, real = real, check = check)


   if(plot == TRUE){
      par(mfrow=c(2,1))
      plot.ts(Re(siginput))
      title("Original signal")
      plot.ts(Re(tmp2$sol))
      title("Reconstructed signal")
   }

   list(sol=tmp2$sol,A=tmp2$A,lam=tmp2$lam,dualwave=tmp2$dualwave,
	gaborets=tmp2$gaborets, solskel=tmp2$solskel,
	inputskel = tmp2$inputskel, Q2 = Q2)

}




RidgeIrregSampling <- function(phinode, node, lam, prob)
#########################################################################
#    RidgeSamplingData:
#    -----------------
#      Given a ridge phi (for the Gabor transform), and lam returns a 
#        (data-driven) subsampled version of length nbnodes.
#
#      Input:
#      ------
#      phinode: ridge at node
#      lam : vector of 2 * number of regular sampled node
#      prob : the percentage of lam to be preserved
#
#      Output:
#      -------
#      node: time coordinates of the ridge sampled by lambda
#      phinode: frequency coordinates of the ridge sampled by lambda
#
#########################################################################
{
   nbnodes <- length(node)   
   mlam <- numeric(nbnodes)
   for(j in 1:nbnodes) 
     mlam[j] <- Mod(lam[j] + 1i*lam[nbnodes + j])

   pct <- quantile(mlam,prob)

   count <- sum(mlam > pct)

   newnode <- numeric(count)
   newphinode <- numeric(count)

   count <- 0
   for(j in 1:nbnodes) {
     if(mlam[j] > pct) {
        newnode[count+1] <- node[j]
        newphinode[count+1] <- phinode[j]             
        count <- count + 1
     }
   }

   list(node = newnode,phinode = newphinode, nbnodes = count)
}





