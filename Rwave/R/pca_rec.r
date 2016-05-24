#########################################################################
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################


#########################################################################
#
#	Functions to reconstruct a signal from the output of
#	 the pca climber algorithm.
#
#########################################################################




pcarec <- function(siginput,inputwt, beemap, orientmap, noct, nvoice, compr,
	maxchnlng=as.numeric(dim(beemap)[1])+10,minnbnodes = 2, 
	w0 = 2*pi, nbchain=100,bstep = 1,ptile =.01,para=5,plot=2,check=FALSE)
#########################################################################
#     pcarec:
#     -------
#      Reconstruction of a real valued signal from ridges found by 
#      pca climbers on a wavelet transform.
#
#      input:
#      ------
#      siginput: input signal
#      inputwt: continuous wavelet transform (output of cwt)
#      beemap: output of pca climber algorithm
#      orientmap: principle direction of pca climber
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      compr: subsampling rate for the ridge
#      maxchnlng: maxlength of a chain of ridges (a,b)
#      bstep: used for the chaining
#      ptile: 
#      para:
#      plot: plot the original and reconstructed signal in display
#      check: check whether the reconstructed signal keeps the values
#             at ridges
#
#      output:
#      -------
#      rec: reconstructed signal
#      ordered: image of the ridges (with different colors)
#      comp: 2D array containing the signals reconstructed from ridges
#
#########################################################################
{

   tmp <- pcafamily(beemap,orientmap,maxchnlng=maxchnlng,bstep=bstep,nbchain=nbchain,ptile=ptile)
   chain <- tmp$chain
   nbchain <- tmp$nbchain
   ordered <- tmp$ordered
   sigsize <- length(siginput)
   rec <- numeric(sigsize)
   plnb <- 0

   if(plot != FALSE){
     par(mfrow=c(2,1))
     plot.ts(siginput)
     title("Original signal")
     image(tmp$ordered)
     title("Chained Ridges")
   }

   sol <- matrix(0,nbchain,sigsize)

   totnbnodes <- 0
   idx <- numeric(nbchain)	
   p <- 0

   if(check==TRUE) {
      inputskel <- matrix(0+0i,nbchain,sigsize)
      solskel <- matrix(0+0i, nbchain, sigsize)
   }

   for (j in 1:nbchain){
      
      sol[j,] <- 0
      nbnode <- chain[j,1]
      bnode <- numeric(nbnode)
      anode <- numeric(nbnode)
      for(k in 1:nbnode) {
        anode[k] <- chain[j,2*k]
        bnode[k] <- chain[j,2*k+1]
      }
      cat("Chain number",j,"\n")

      tmp2 <- pcaregrec(siginput[min(bnode):max(bnode)],
          inputwt[min(bnode):max(bnode),],
          anode,bnode,compr,noct,nvoice, 
          w0 = w0, para = para,minnbnodes = minnbnodes, check=check);

      if(is.list(tmp2)==TRUE) {
        totnbnodes <- totnbnodes + tmp2$nbnodes
        if((sigsize-min(bnode)) > (length(tmp2$sol)-tmp2$bstart)){
          np <- length(tmp2$sol) - tmp2$bstart
          sol[j,min(bnode):(np+min(bnode))]<-tmp2$sol[tmp2$bstart:length(tmp2$sol)]
        }
        else {
	  np <- sigsize - min(bnode)
          sol[j,min(bnode):sigsize]<-tmp2$sol[tmp2$bstart:(np+tmp2$bstart)]
        }
        if(min(bnode) < tmp2$bstart) {
          np <- min(bnode)-1
          sol[j,1:min(bnode)]<-tmp2$sol[(tmp2$bstart-np):(tmp2$bstart)]
        }          
        else {
          np <- tmp2$bstart-1
          sol[j,(min(bnode)-np):min(bnode)]<-tmp2$sol[1:(tmp2$bstart)]
        }

        if(check==TRUE) {
           bridge <- tmp2$bnode
           aridge <- tmp2$anode
           wtsol <- cwt(sol[j,],noct,nvoice)
           for(k in 1:length(bridge)) solskel[j,k] <- wtsol[bridge[k],aridge[k]]
           for(k in 1:length(bridge)) inputskel[j,k] <- inputwt[bridge[k],aridge[k]]
        }
        rec <- rec+sol[j,]
      }
      plnb <- plnb + 1
      p <- p + 1
      idx[p] <- j
   }

   if(plot == 1){
      par(mfrow=c(2,1))
      par(cex=1.1)
      plot.ts(siginput)
      title("Original signal")
      plot.ts(Re(rec))
	
      title("Reconstructed signal")
   }
   else if (plot == 2){
      par(mfrow=c(plnb+2,1))
      par(mar=c(2,4,4,4))
      par(cex=1.1)
      par(err=-1)
      plot.ts(siginput)
      title("Original signal")

      for (j in 1:p)
        plot.ts(sol[idx[j],]);
      plot.ts(Re(rec))	
      title("Reconstructed signal")
   }

   cat("Total number of ridge samples used: ",totnbnodes,"\n")

   par(mfrow=c(1,1))
   if(check==TRUE) list(rec=rec,ordered=ordered,chain=chain,
                     comp=tmp,inputskel=inputskel,solskel=solskel,lam=tmp2$lam)
   else list(rec=rec, ordered=ordered,chain=chain,comp=tmp)
}

PcaRidgeSampling <- function(anode, bnode, compr)
#########################################################################
#    PcaRidgeSampling:
#    ----------------
#      Given a ridge (a,b)(for the Gabor transform), returns a 
#        (regularly) subsampled version of length nbnodes.
#
#      Input:
#      ------
#      anode: scale coordinates of the (unsampled) ridge
#      bnode: time coordinates of the (unsampled) ridge
#      compr: compression ratio to obtain from the ridge
#
#      Output:
#      -------
#      bnode: time coordinates of the ridge samples (from 1 to nbnodes
#             step by compr)
#      anode: scale coordinates of the ridge samples
#      nbnode: number of node sampled
#
#########################################################################
{
   # number of node at the ridge
   nbnode <- length(anode)

   # compr to be integer
   compr <- as.integer(compr)
   if(compr < 1) compr <- 1

   # sample rate is compr at the ridge
   k <- 1
   count <- 1
   while(k < nbnode) {
     anode[count] <- anode[k]
     bnode[count] <- bnode[k]
     k <- k + compr 
     count <- count + 1
   }

   anode[count] <- anode[nbnode]
   bnode[count] <- bnode[nbnode]
   
   aridge <- numeric(count)
   bridge <- numeric(count)
   aridge <- anode[1:count]
   bridge <- bnode[1:count]
   
   list(anode = aridge, bnode = bridge, nbnodes = count)
}


pcaregrec <- function(siginput,cwtinput,anode,bnode,compr,noct,nvoice,
	w0 = 2*pi, plot = FALSE, para = 5, minnbnodes = 2, check=FALSE)
#########################################################################
#     pcaregrec:
#     ---------
#      Reconstruction of a real valued signal from a (continuous) ridge
#      (uses a regular sampling of the ridge)
#
#      input:
#      -------
#      siginput: input signal
#      cwtinput: Continuous wavelet transform (output of cwt)
#      anode: the scale coordinates of the (unsampled) ridge
#      bnode: the time coordinates of the (unsampled) ridge
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
#      Q2: second part of the reconstruction kernel
#      nbnodes: number of nodes used for the reconstruction.
#
##########################################################################
{

#  Generate Sampled ridge
#   tmp <- wRidgeSampling(anode,compr,nvoice)
#   bnode <- tmp$node
#   anode <- tmp$phinode
#   nbnodes <- tmp$nbnodes

   tmp <- PcaRidgeSampling(anode, bnode, compr)
   bnode <- tmp$bnode
   anode <- tmp$anode
   nbnodes <- tmp$nbnodes

   cat("Sampled nodes (b,a): \n")
   for(j in 1:nbnodes) cat("(",bnode[j],anode[j],")")
   cat("\n")

   if(nbnodes < minnbnodes){
      cat(" Chain too small\n")
      NULL
   }
   else {

    a.min <- 2 * 2^(min(anode)/nvoice)
    a.max <- 2 * 2^(max(anode)/nvoice)
    b.min <- min(bnode)
    b.max <- max(bnode)
    b.max <- (b.max-b.min+1) + round(para * a.max)
    b.min <- (b.min-b.min+1) - round(para * a.max)

# location 2-b.min is the place where the min(bnode) is.
# The position at 2-b.min of returned signal is aligned to 
# the position at min(b.node) for the reconstruction

    bnode <- bnode-min(bnode)+2-b.min
    b.inc <- 1
    np <- as.integer((b.max - b.min)/b.inc) +1
    cat("(size:",np,",",nbnodes,"sampled nodes):\n")

    # Generating the Q2 term in reconstruction kernel
    Q2 <- 0
	
    # Generating the Q1 term in reconstruction kernel
    one <- numeric(np)
    one[] <- 1
    Qinv <-  1/one

    tmp2 <- pcaridrec(cwtinput,bnode,anode,noct,nvoice,
              Qinv,np,w0=w0,check=check)

    if(plot == TRUE){
       par(mfrow=c(2,1))
       plot.ts(Re(siginput))
       title("Original signal")
       plot.ts(Re(tmp2$sol))
       title("Reconstructed signal")
    }
    lam <- tmp2$lam
#    if(check==TRUE) {
#      N <- length(lam)/2
#      mlam <- numeric(N)
#      for(j in 1:N)
#         mlam[j] <- Mod(lam[j] + i * lam[N + j])
#      npl(2)
#      plot.ts(lam,xlab="Number",ylab="Lambda Value")
#      title("Lambda Profile")
#      plot.ts(sort(mlam),xlab="Number",ylab="Mod of Lambda")
#    }
    list(sol = tmp2$sol,A = tmp2$A,
          lam = tmp2$lam, dualwave = tmp2$dualwave, 
          morvelets = tmp2$morvelets,pQ2 = Q2, nbnodes=nbnodes,
          bstart=2-b.min, anode=tmp$anode,bnode=tmp$bnode)
   }
	
}


pcaridrec <- function(cwtinput,bnode,anode,noct,nvoice,Qinv,np,
                      w0 = 2*pi, check = FALSE)
#########################################################################
#     pcaridrec:
#     ---------
#      Reconstruction of a real valued signal from a ridge, given
#       the kernel of the bilinear form.
#
#      input:
#      ------
#      cwtinput: Continuous wavelet transform (output of cwt)
#      bnode: time coordinates of the ridge samples
#      anode: scale coordinates of the ridge samples
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      Qinv: inverse of the reconstruction kernel
#      epsilon: coefficient of the Q2 term in the reconstruction kernel
#      np: number of samples of the reconstructed signal
#      w0: central frequency of Morlet wavelet
#      check: if set to TRUE, computes the wavelet transform of the
#            reconstructed signal
#
#      output:
#      -------
#      sol: reconstruction from a ridge
#      A: <wavelets,dualwavelets> matrix 
#      lam: coefficients of dual wavelets in reconstructed signal.
#      morwave: array of morlet wavelets located on the ridge samples
#      dualwave: array containing the dual wavelets.
#      solskel: wavelet transform of sol, restricted to the ridge
#      inputskel: wavelet transform of signal, restricted to the ridge
#
#########################################################################
{
   # number of sampled nodes in a ridge
   N <- length(bnode)

   aridge <- anode
   bridge <- bnode

   morvelets <- pcamorwave(bridge,aridge,nvoice,np,N)

   cat("morvelets;  ")

   sk <- pcazeroskeleton(cwtinput,Qinv,morvelets,bridge,aridge,N)    

   cat("skeleton.\n")

   solskel <- 0 #not needed if check not done
   inputskel <- 0 #not needed if check not done

   list(sol=sk$sol,A=sk$A,lam=sk$lam,dualwave=sk$dualwave,morvelets=morvelets,
	solskel=solskel,inputskel = inputskel)
}


pcazeroskeleton <- function(cwtinput,Qinv,morvelets,bridge,aridge,N)
#########################################################################
#     pcazeroskeleton:
#     ----------------
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
   constraints2 <- dim(morvelets)[2]
   dualwave <- matrix(0,tmp1,constraints2)
   for (j in 1:tmp1) {
      dualwave[j,] <- morvelets[j,]*Qinv[j]
   }
   A <- t(morvelets) %*% dualwave

   rskel <- numeric(2*N)

   if(is.vector(cwtinput) == TRUE) {
     for(j in 1:N) {
      rskel[j] <- Re(cwtinput[aridge[j]])
      rskel[N+j] <- -Im(cwtinput[aridge[j]])
     }
   }
   else {
     for(j in 1:N) {
      rskel[j] <- Re(cwtinput[bridge[j]-bridge[1]+1,aridge[j]])
      rskel[N+j] <- -Im(cwtinput[bridge[j]-bridge[1]+1,aridge[j]])
     }
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


pcamorwave <- function(bridge, aridge, nvoice, np, N, w0 = 2*pi)
#########################################################################
#     pcamorwave: ( the same function of morwave )
#     -----------
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
#      cwtdirac <- vwt(dirac,scale);
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



