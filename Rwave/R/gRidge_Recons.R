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
#                  Princeton University 
#                  All right reserved                           
#########################################################################





gregrec <- function(siginput, gtinput, phi, nbnodes, nvoice,
	freqstep, scale, epsilon = 0, fast = FALSE, plot = FALSE,
	para = 0, hflag = FALSE, real = FALSE, check = FALSE)
#########################################################################
#     gregrec:
#     -------
#      Reconstruction of a real valued signal from a (continuous) 
#      Gabor ridge (uses a regular sampling of the ridge)
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
#      para: scale parameter for extrapolating the ridges.
#      hflag:if set to TRUE, uses $Q_1$ as first term in the kernel.
#      real: if set to TRUE, only uses constraints on the real part
#            of the transform for the reconstruction.
#      check: if set to TRUE, computes the wavelet transform of the
#            reconstructed signal
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <gaborlets,dualgaborlets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#      gaborets: array containing the wavelets on sampled ridge.
#      solskel: Gabor transform of sol, restricted to the ridge
#      inputskel: Gabor transform of signal, restricted to the ridge
#      Q2: second part of the reconstruction kernel
#########################################################################
{
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

   if(plot == TRUE){
      par(mfrow=c(2,1))
      plot.ts(Re(siginput))
      title("Original signal")
      plot.ts(Re(tmp2$sol))
      title("Reconstructed signal")
   }

npl(2)
lam <- tmp2$lam
plot.ts(lam,xlab="Number",ylab="Lambda Value")
title("Lambda Profile")
N <- length(lam)/2
mlam <- numeric(N)
for(j in 1:N)
   mlam[j] <- Mod(lam[j] + 1i*lam[N + j])
plot.ts(sort(mlam))

   list(sol=tmp2$sol,A=tmp2$A,lam=tmp2$lam,dualwave=tmp2$dualwave,
	gaborets=tmp2$gaborets, solskel=tmp2$solskel,
	inputskel = tmp2$inputskel, Q2 = Q2)

}




gridrec <- function(gtinput, node, phinode, nvoice,
	freqstep, scale, Qinv, epsilon, np, real = FALSE, check = FALSE)
#########################################################################
#     gridrec:
#     --------
#      Reconstruction of a real valued signal from a gabor ridge.
#
#      input:
#      ------
#      gtinput: Continuous gabor transform (output of cgt)
#      node: time coordinates of the ridge samples.
#      phinode: frequency coordinates of the ridge samples.
#      nvoice: number of different frequencies.
#      freqstep: sampling rate for the frequency axis.
#      scale: scale of the window.
#      Qinv: inverse of the reconstruction kernel
#      epsilon: coefficient of the Q2 term in the reconstruction kernel
#      np: number of samples of the reconstructed signal
#      real: if set to TRUE, only uses constraints on the real part
#            of the transform for the reconstruction.
#      check: if set to TRUE, computes the Gabor transform of the
#            reconstructed signal.
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <gaborlets,dualgaborlets> matrix.
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual gaborlets.
#      gaborets: array of gaborlets located on the ridge samples
#      solskel: Gabor transform of sol, restricted to the ridge
#      inputskel: Gabor transform of signal, restricted to the ridge
#
#########################################################################
{
   N <- length(node)


   omegaridge <- phinode
   bridge <- node
   if(real == TRUE)
      gaborets <- gwave2(bridge,omegaridge,nvoice,freqstep,scale,np,N)
   else
      gaborets <- gwave(bridge,omegaridge,nvoice,freqstep,scale,np,N)

   cat("gaborets;")

   if(epsilon == 0){
      if(real == TRUE)
         sk <- zeroskeleton2(gtinput,Qinv,gaborets,bridge,omegaridge,N)
      else
         sk <- zeroskeleton(gtinput,Qinv,gaborets,bridge,omegaridge,N)
      }

   else
      sk <- skeleton(gtinput,Qinv,gaborets,bridge,omegaridge,N)

   cat("skeleton.\n")

   solskel <- 0 #not needed if check not done
   inputskel <- 0 #not needed if check not done

   if(check == TRUE){
      solskel <- complex(N)
      inputskel <- complex(N)
      gtsol <- cgt(sk$sol,nvoice,freqstep,plot=FALSE)
      for(j in 1:N) solskel[j] <- gtsol[bridge[j],omegaridge[j]]

      for (j in 1:N)
         inputskel[j] <- gtinput[bridge[j],omegaridge[j]]
   }

   list(sol=sk$sol,A=sk$A,lam=sk$lam,dualwave=sk$dualwave,gaborets=gaborets,
	solskel=solskel,inputskel = inputskel)
}


gwave <- function(bridge,omegaridge,nvoice,freqstep,scale,np,N)
#########################################################################
#     gwave:
#     ------
#      Generation of the gaborets located on the ridge
#
#      input:
#      ------
#      bridge: time coordinates of the ridge samples
#      omegaridge: frequency coordinates of the ridge samples
#      nvoice: number of different scales per octave
#      freqstep: sampling rate for the frequency axis
#      scale: scale of the window
#      np: size of the reconstruction kernel
#      N: number of complex constraints
#
#      output:
#      -------
#      gaborets: array of morlet wavelets located on the ridge samples
#
#########################################################################
{

   gaborets <- matrix(0,np,2*N)

    omegaridge <- omegaridge * freqstep
    tmp <- vecgabor(np,N,bridge,omegaridge,scale)
    dim(tmp) <- c(np,N)
    gaborets[,1:N] <- Re(tmp[,1:N])
    gaborets[,(N+1):(2*N)] <- Im(tmp[,1:N])

# Alternative way of generating gaborets
#   for(k in 1:N)  {
#      omega <- omegaridge[k]*freqstep;
#      tmp <- gabor(np,bridge[k],omega,scale);
#      gaborets[,k] <- Re(tmp);
#      gaborets[,k+N] <- Im(tmp);
#    }


# Still another way of generating gaborets
#
#   dirac <- 1:np
#   dirac[] <- 0
#   for(k in 1:N)  {
#      dirac[bridge[k]] <- 1.0;
#      omega <- omegaridge[k]*freqstep;
#      gtdirac <- vgt(dirac,omega,scale,open=FALSE,plot=FALSE);
#      gaborets[,k] <- Re(gtdirac);
#      gaborets[,k+N] <- Im(gtdirac);
#      dirac[] <- 0
#    }

    gaborets
}



gwave2 <- function(bridge,omegaridge,nvoice,freqstep,scale,np,N)
#########################################################################
#     gwave2:
#     -------
#      Generation of the real parts of gaborets located on the ridge
#
#      input:
#      ------
#      bridge: time coordinates of the ridge samples
#      omegaridge: frequency coordinates of the ridge samples
#      nvoice: number of different scales per octave
#      freqstep: sampling rate for the frequency axis
#      scale: scale of the window
#      np: size of the reconstruction kernel
#      N: number of complex constraints
#
#      output:
#      -------
#      gaborets: array of morlet wavelets located on the ridge samples
#
#########################################################################
{
    omegaridge <- omegaridge * freqstep
    gaborets <- Re(vecgabor(np,N,bridge,omegaridge,scale))
    dim(gaborets) <- c(np,N)


#   for(k in 1:N)  {
#      omega <- omegaridge[k]*freqstep;
#      tmp <- gabor(np,bridge[k],omega,scale);
#      gaborets[,k] <- Re(tmp);
#    }
    gaborets
}




gsampleOne <- function(node,scale,np)
#########################################################################
#     gsampleOne:
#     -----------
#
#       Generate a ``sampled identity"
#
#     input:
#     ------
#      node: location of the reconstruction gabor functions
#      scale: scale of the gabor functions
#      np: size of the reconstructed signal
#
#
#     output:
#     -------
#      dia: diagonal of the ``sampled'' Q1 term (1D vector)
#
#########################################################################
{
   dia <- numeric(np)
   N <- length(node)

   normalization <- sqrt(pi) * scale

   for(j in 1:np) {
     tmp1 <- (j-node)/scale
     tmp1 <- exp(-(tmp1 * tmp1))

     dia[j] <- sum(tmp1) #/normalization
    }
   dia

#   for(j in 1:np) {
#     tmp <- 0
#     for(k in 1:N) {
#       b <- node[k]
#       tmp1 <- (j-b)/scale
#       tmp1 <- exp(-(tmp1 * tmp1))
#       tmp <- tmp + tmp1
#     }
#     dia[j] <- tmp #/normalization
#   }
   dia
}






