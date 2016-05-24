#########################################################################
#      $Log: Ridge_Recons.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by     
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





regrec <- function(siginput, cwtinput, phi, compr, noct, nvoice,
                   epsilon = 0, w0 = 2*pi, fast = FALSE, plot = FALSE,
                   para = 0, hflag = FALSE, check = FALSE, minnbnodes = 2,
                   real = FALSE)
#########################################################################
#     regrec:
#     -------
#      Reconstruction of a real valued signal from a (continuous) ridge
#      (uses a regular sampling of the ridge)
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
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#      morvelets: array containing the wavelets on sampled ridge.
#      solskel: wavelet transform of sol, restricted to the ridge
#      inputskel: wavelet transform of signal, restricted to the ridge
#      Q2: second part of the reconstruction kernel
#      nbnodes: number of nodes used for the reconstruction.
#
##########################################################################
{
   ## Generate Sampled ridge

   tmp <- wRidgeSampling(phi,compr,nvoice)

   node <- tmp$node
   phinode <- tmp$phinode
   nbnodes <- tmp$nbnodes
   phinode <- as.integer(phinode)	

   if(nbnodes < minnbnodes){
     cat(" Chain too small\n")
     NULL
   }
   else {
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

   fast <- FALSE

   ## Generating the Q2 term in reconstruction kernel
   if(epsilon == 0)
     Q2 <- 0
   else {
     if (fast == FALSE)
       Q2 <- rkernel(node, phinode, nvoice, x.min = x.min, x.max = x.max,
                     w0 = w0)
     else
       Q2 <- fastkernel(node, phinode, nvoice, x.min = x.min,
                        x.max = x.max, w0 = w0)
   }
   
   ## Generating the Q1 term in reconstruction kernel
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
   
   tmp2 <- ridrec(cwtinput, node, phinode, noct, nvoice,
                  Qinv, epsilon, np, w0=w0, check=check, real=real)
   
   if(plot == TRUE) {
     par(mfrow=c(2,1))
     plot.ts(Re(siginput))
     title("Original signal")
     plot.ts(Re(tmp2$sol))
     title("Reconstructed signal")
   }
   npl(2)
   lam <- tmp2$lam
   if(plot == TRUE) {
     plot.ts(lam,xlab="Number",ylab="Lambda Value")
     title("Lambda Profile")
   }
   N <- length(lam)/2
   mlam <- numeric(N)
   for(j in 1:N)
     mlam[j] <- Mod(lam[j] + lam[N + j]*(1i))
   if(plot == TRUE) plot.ts(sort(mlam))
   
   list(sol = tmp2$sol, A = tmp2$A, lam = tmp2$lam, dualwave = tmp2$dualwave,
        morvelets = tmp2$morvelets, solskel = tmp2$solskel,
        inputskel = tmp2$inputskel, Q2 = Q2, nbnodes = nbnodes)
 }
}


ridrec <- function(cwtinput, node, phinode, noct, nvoice,
	Qinv, epsilon, np, w0 = 2*pi, check = FALSE, real = FALSE)
#########################################################################
#     ridrec:
#     ------
#      Reconstruction of a real valued signal from a ridge, given
#       the kernel of the bilinear form.
#
#      input:
#      ------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      node: time coordinates of the ridge samples
#      phinode: scale coordinates of the ridge samples
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      Qinv: inverse of the reconstruction kernel
#      epsilon: coefficient of the Q2 term in the reconstruction kernel
#      np: number of samples of the reconstructed signal
#      w0: central frequency of Morlet wavelet
#      check: if set to TRUE, computes the wavelet transform of the
#            reconstructed signal
#      real: if set to TRUE, only uses constraints on the real part
#            of the transform for the reconstruction.
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#      morvelets: array of morlet wavelets located on the ridge samples
#      solskel: wavelet transform of sol, restricted to the ridge
#      inputskel: wavelet transform of signal, restricted to the ridge
#
#########################################################################
{
   N <- length(node)

   aridge <- phinode
   bridge <- node

   if (real == TRUE)
      morvelets <- morwave2(bridge,aridge,nvoice,np,N)
   else
      morvelets <- morwave(bridge,aridge,nvoice,np,N)

   cat("morvelets;  ")


   if( epsilon == 0){
      if (real == TRUE)
         sk <- zeroskeleton2(cwtinput,Qinv,morvelets,bridge,aridge,N)
      else
         sk <- zeroskeleton(cwtinput,Qinv,morvelets,bridge,aridge,N)      
   }
   else { 
      if(real == TRUE) 
         sk <- skeleton2(cwtinput,Qinv,morvelets,bridge,aridge,N)
      else
         sk <- skeleton(cwtinput,Qinv,morvelets,bridge,aridge,N)
   }
   cat("skeleton.\n")

   solskel <- 0 #not needed if check not done
   inputskel <- 0 #not needed if check not done

   if(check == TRUE){
      wtsol <- cwt(Re(sk$sol),noct,nvoice)
      solskel <- complex(N)
      for(j in 1:N) solskel[j] <- wtsol[bridge[j],aridge[j]]

      inputskel <- complex(N)
      for (j in 1:N)
         inputskel[j] <- cwtinput[bridge[j],aridge[j]]
   }

   list(sol=sk$sol,A=sk$A,lam=sk$lam,dualwave=sk$dualwave,morvelets=morvelets,
	solskel=solskel,inputskel = inputskel)
}



morwave <- function(bridge, aridge, nvoice, np, N, w0 = 2*pi)
#########################################################################
#     morwave:
#     -------
#      Generation of the wavelets located on the ridge.
#
#      input:
#      ------
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      nvoice: number of different scales per octave
#      np: size of the reconstruction kernel
#      N: number of complex constraints
#      w0: central frequency of Morlet wavelet
#
#      output:
#      -------
#      morvelets: array of morlet wavelets located on the ridge samples
#
#########################################################################
{

   morvelets <- matrix(0,np,2*N)

   aridge <-  2 * 2^((aridge - 1)/nvoice)
   tmp <- vecmorlet(np,N,bridge,aridge,w0 = w0)
   dim(tmp) <- c(np,N)
   morvelets[,1:N] <- Re(tmp[,1:N])
   morvelets[,(N+1):(2*N)] <- Im(tmp[,1:N])

#  Another way of generating morvelets
#   dirac <- 1:np
#   dirac[] <- 0
#   for(k in 1:N)  {
#      dirac[bridge[k]] <- 1.0;
#      scale <- 2 * 2^((aridge[k]-1)/nvoice);
#      cwtdirac <- vwt(dirac,scale,open=FALSE);
#      morvelets[,k] <- Re(cwtdirac);
#      morvelets[,k+N] <- Im(cwtdirac);
#      dirac[] <- 0
#    }

#  Still another way of generating morvelets
#   for(k in 1:N)  {
#      scale <- 2 * 2^((aridge[k]-1)/nvoice);
#      tmp <- morlet(np,bridge[k],scale)/sqrt(2*pi);
#      morvelets[,k] <- Re(tmp);
#      morvelets[,k+N] <- Im(tmp);
#    }

    morvelets
}



morwave2 <- function(bridge, aridge, nvoice, np, N, w0 = 2*pi)
#########################################################################
#     morwave2:
#     ---------
#      Generation of the real parts of morlet wavelets
#        located on the ridge.
#
#      input:
#      -------
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      nvoice: number of different scales per octave
#      np: size of the reconstruction kernel
#      N: number of real constraints
#      w0: central frequency of Morlet wavelet
#
#      output:
#      -------
#      morvelets: array of morlet wavelets located on the ridge samples
#
#########################################################################
{

   morvelets <- matrix(0,np,N)

   aridge <-  2 * 2^((aridge - 1)/nvoice)
   morvelets <- Re(vecmorlet(np,N,bridge,aridge,w0 = w0))
   dim(morvelets)<- c(np,N)

# changed by Wen 5/22
#
#   aridge <-  2 * 2^((aridge - 1)/nvoice)
#   morvelets <- Re(vecmorlet(np,N,bridge,aridge,w0 = w0))
#   dim(morvelets) <- c(np,N)

   morvelets
}



skeleton <- function(cwtinput, Qinv, morvelets, bridge, aridge, N)
#########################################################################
#     skeleton:
#     --------
#      Computation of the reconstructed signal from the ridge
#
#      input:
#      -------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      Qinv: inverse of the reconstruction kernel (2D array)
#      morvelets: array of morlet wavelets located on the ridge samples
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      N: number of complex constraints
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#
##########################################################################
{
   dualwave <- Qinv %*% morvelets
   A <- t(morvelets) %*% dualwave

   rskel <- numeric(2*N)
   for(j in 1:N) {
      rskel[j] <- Re(cwtinput[bridge[j]-bridge[1]+1,aridge[j]]);
      rskel[N+j] <- -Im(cwtinput[bridge[j]-bridge[1]+1,aridge[j]])
   }

  B <- SVD(A)
  d <- B$d
  Invd <- numeric(length(d))
  for(j in 1:(2*N))
    if(d[j] < 1.0e-6) 
       Invd[j] <- 0
    else Invd[j] <- 1/d[j]

  lam <- B$v %*% diag(Invd) %*% t(B$u) %*% rskel
   sol <- dualwave%*%lam

   list(lam=lam,sol=sol,dualwave=dualwave,A=A)
}


skeleton2 <- function(cwtinput, Qinv, morvelets, bridge, aridge, N)
#########################################################################
#     skeleton2:
#     ---------
#      Computation of the reconstructed signal from the ridge, in the
#        case of real constraints.
#
#      input:
#      -------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      Qinv: inverse of the reconstruction kernel (2D array)
#      morvelets: array of morlet wavelets located on the ridge samples
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      N: number of real constraints
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#
##########################################################################
{
   dualwave <- Qinv %*% morvelets
   A <- t(morvelets) %*% dualwave

   rskel <- numeric(N)
   for(j in 1:N) {
      rskel[j] <- Re(cwtinput[bridge[j]-bridge[1]+1,aridge[j]]);
   }
   B <- SVD(A)
   d <- B$d
   Invd <- numeric(length(d))
   for(j in 1:N)
     if(d[j] < 1.0e-6) 
       Invd[j] <- 0
    else Invd[j] <- 1/d[j]

   lam <- B$v %*% diag(Invd) %*% t(B$u) %*% rskel

   sol <- dualwave%*%lam

   list(lam=lam,sol=sol,dualwave=dualwave,A=A)
}



zeroskeleton <- function(cwtinput, Qinv, morvelets, bridge, aridge, N)
#########################################################################
#     zeroskeleton:
#     -------------
#      Computation of the reconstructed signal from the ridge when
#       the epsilon parameter is set to 0.
#
#      input:
#      -------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      Qinv: 1D array (warning: different from skeleton) containing
#            the diagonal of the inverse of reconstruction kernel.
#      morvelets: array of morlet wavelets located on the ridge samples
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      N: number of complex constraints
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#
##########################################################################
{
   tmp1 <- dim(morvelets)[1]
   tmp2 <- dim(morvelets)[2]
cat(dim(morvelets))
   dualwave <- matrix(0,tmp1,tmp2)
   for (j in 1:tmp1)
      dualwave[j,] <- morvelets[j,]*Qinv[j]
   A <- t(morvelets) %*% dualwave

   rskel <- numeric(2*N)

   for(j in 1:N) {
      rskel[j] <- Re(cwtinput[bridge[j]-bridge[1]+1,aridge[j]]);
      rskel[N+j] <- -Im(cwtinput[bridge[j]-bridge[1]+1,aridge[j]])
   }

   B <- SVD(A)
   d <- B$d
   Invd <- numeric(length(d))
   for(j in 1:(2*N))
     if(d[j] < 1.0e-6) 
       Invd[j] <- 0
    else Invd[j] <- 1/d[j]

   lam <- B$v %*% diag(Invd) %*% t(B$u) %*% rskel

   sol <- dualwave%*%lam

   list(lam=lam,sol=sol,dualwave=dualwave,A=A)
}




zeroskeleton2 <- function(cwtinput, Qinv, morvelets, bridge, aridge, N)
#########################################################################
#     zeroskeleton2:
#     -------------
#      Computation of the reconstructed signal from the ridge when
#       the epsilon parameter is set to 0 and the constraints are real.
#
#      input:
#      -------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      Qinv: 1D array (warning: different than in skeleton) containing
#            the diagonal of the inverse of reconstruction kernel.
#      morvelets: array of morlet wavelets located on the ridge samples
#      bridge: time coordinates of the ridge samples
#      aridge: scale coordinates of the ridge samples
#      N: number of real constraints
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#
##########################################################################
{
   tmp1 <- dim(morvelets)[1]
   tmp2 <- dim(morvelets)[2]
   dualwave <- matrix(0,tmp1,tmp2)
   for (j in 1:tmp1)
      dualwave[j,] <- morvelets[j,]*Qinv[j]
   A <- t(morvelets) %*% dualwave

   rskel <- numeric(N)

   for(j in 1:N) 
     rskel[j] <- Re(cwtinput[bridge[j]-bridge[1]+1,aridge[j]]);

   B <- SVD(A)
   d <- B$d
   Invd <- numeric(length(d))
   for(j in 1:N)
     if(d[j] < 1.0e-6) 
       Invd[j] <- 0
    else Invd[j] <- 1/d[j]
   lam <- B$v %*% diag(Invd) %*% t(B$u) %*% rskel

   sol <- dualwave%*%lam

   list(lam=lam,sol=sol,dualwave=morvelets,A=A)
}





regrec2 <- function(siginput, cwtinput, phi, nbnodes, noct, nvoice, Q2,
	epsilon = 0.5, w0 = 2*pi , plot = FALSE)
#########################################################################
#     regrec2:
#     -------
#      Reconstruction of a real valued signal from a (continuous) ridge
#      (uses a regular sampling of the ridge), from a precomputed kernel.
#
#      input:
#      -------
#      siginput: input signal
#      cwtinput: Continuous wavelet transform (output of cwt)
#      phi: (unsampled) ridge
#      nbnodes: number of samples on the ridge
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      Q2: second part of the reconstruction kernel
#      epsilon: coeff of the Q2 term in reconstruction kernel
#      w0: central frequency of Morlet wavelet
#      plot: if set to TRUE, displays original and reconstructed signals
#
#      output: same as ridrec
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      dualwave: array containing the dual wavelets.
#      morvelets: array of morlet wavelets located on the ridge samples
#      solskel: wavelet transform of sol, restricted to the ridge
#      inputskel: wavelet transform of signal, restricted to the ridge
#
#########################################################################
{
   tmp <- RidgeSampling(phi,nbnodes)
   node <- tmp$node
   phinode <- tmp$phinode

   np <- dim(Q2)[1]
   one <- diag(np)
   Q <- one + epsilon * Re(Q2)

   if(epsilon == 0)
      Qinv <- one
   else
      Qinv <- solve(Q)

   tmp2 <- ridrec(cwtinput,node,phinode,noct,nvoice,Qinv,epsilon,np)

   if(plot == TRUE){
      par(mfrow=c(2,1))
      plot.ts(Re(siginput))
      title("Original signal")
      plot.ts(Re(tmp2$sol))
      title("Reconstructed signal")
   }
   tmp2
}




RidgeSampling <- function(phi, nbnodes)
#########################################################################
#    RidgeSampling:
#    --------------
#      Given a ridge phi (for the Gabor transform), returns a 
#        (regularly) subsampled version of length nbnodes.
#
#      Input:
#      ------
#      phi: ridge
#      nbnodes: number of samples.
#
#      Output:
#      -------
#      node: time coordinates of the ridge samples
#      phinode: frequency coordinates of the ridge samples
#
#########################################################################
{
   node <- numeric(nbnodes)
   phinode <- numeric(nbnodes)
   N <- nbnodes
   LL <- length(phi)

   for (i in 1:(N+1)) node[i] <- as.integer(((LL-1)*i+nbnodes-LL+1)/nbnodes)
   for (i in 1:(N+1)) phinode[i] <- phi[node[i]]

   list(node = node,phinode = phinode)
}



wRidgeSampling <- function(phi, compr, nvoice)
#########################################################################
#    wRidgeSampling:
#    ---------------
#      Given a ridge phi, returns a subsampled version of length
#      nbnodes (wavelet case).
#
#      Input:
#      ------
#      phi: ridge
#      compr: subsampling rate for the wavelet coefficients (at scale 1)
#             when compr <= 0, uses all the points in a ridge
#
#      Output:
#      -------
#      node: time coordinates of the ridge samples
#      phinode: frequency coordinates of the ridge samples
#      nbnodes: number of samples.
#
#########################################################################
{


   LL <- length(phi)
   node <- numeric(LL)
   phinode <- numeric(LL)

   LN <- 1
   j <- 1
   node[1] <- 1

# We multiply a number on the phi. This number is obtained experimentally
# If we have a ridge which is linear of scale, then the subsampling is
# adapted with s, however, in the wavelet case, the ridge is linear 
# in log(scale), the subsampling rate is therefore adapted to k*s^l
# l here is the constant 100/15, and k is compr 
   aridge <- 2 * 2^((phi-1)/(nvoice * 100/15))

#   while (node[j] < LL){
#      j <- j + 1
#      LN <- round(compr * aridge[node[j-1]])
#      cat("LN=",LN)
#      node[j] <- node[j-1] + LN
#      cat("; node=",node[j],"\n")
#      phinode[j] <- phi[node[j]]
#      cat("; phinode=",phinode[j],"\n")
#   }

   nbnodes <- 0
   while (node[j] < LL){
      phinode[j] <- phi[node[j]]
      nbnodes <- nbnodes + 1

# To guarantee at least advance by one 
       LN <- max(1,round(compr * aridge[node[j]]))

      j <- j + 1
      node[j] <- node[j-1] + LN
   }


   tmpnode <- numeric(nbnodes)
   tmpphinode <- numeric(nbnodes)
   tmpnode <- node[1:nbnodes]
   tmpphinode <- phinode[1:nbnodes]

   list(node = tmpnode,phinode = tmpphinode, nbnodes = nbnodes)
}




sridrec <- function(tfinput, ridge)
#########################################################################
#     sridrec
#     -------
#      Simple reconstruction of a real valued signal from a ridge,
#      by restriction of the wavelet transform.
#
#      Input:
#      ------
#     tfinput: Continuous time-frequency transform (output of cwt or cgt)
#     ridge: (unsampled) ridge
#
#      Output:
#      -------
#      rec: reconstructed signal
#
#########################################################################
{
   bsize <- dim(tfinput)[1]
   asize <- dim(tfinput)[2]
  
   rec <- numeric(bsize)
   rec[] <- 0 + 0i

   for(i in 1:bsize) {
      if(ridge[i] > 0)
        rec[i] <- rec[i] + 2 * tfinput[i,ridge[i]]
   }
  
   Re(rec)
}




