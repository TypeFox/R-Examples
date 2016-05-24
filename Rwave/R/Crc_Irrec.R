#########################################################################
#      $Log: Crc_Rec.S,v $
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


#########################################################################
#
#	Functions to reconstruct a signal from the output of
#	 the crazy climber algorithm.
#
#########################################################################

crcirrec <- function(siginput, inputwt, beemap, noct, nvoice, compr,
                     minnbnodes = 2, w0 = 2*pi, bstep = 5, ptile = .01,
                     prob = 0.8, epsilon = .5, fast = FALSE, para = 0,
                     real = FALSE, plot=1)
#########################################################################
#     crcirrec:
#     -------
#      Reconstruction of a real valued signal from ridges found by 
#      crazy climbers on a wavelet transform (sampled irregularly).
#
#      input:
#      ------
#      siginput: input signal
#      inputwt: continuous wavelet transform (output of cwt)
#      beemap: output of crazy climber algorithm
#      noct: number of octaves (powers of 2)
#      nvoice: number of different scales per octave
#      compr: subsampling rate for the ridge
#      bstep: used for the chaining
#      ptile: 
#      epsilon: coeff of the Q2 part of reconstruction kernel
#      fast: if set to TRUE, computes the Q2 kernel by Riemann sums
#            if not, uses a Romberg adaptive step quadrature.
#      para:
#      plot: plot the original and reconstructed signal in display
#      prob: the percentile to keep 
#
#      output:
#      -------
#      rec: reconstructed signal
#      ordered: image of the ridges (with different colors)
#      comp: 2D array containing the signals reconstructed from ridges
#
#########################################################################
{
  tmp <- cfamily(beemap,bstep,ptile=ptile)
  image(tmp$ordered)
  chain <- tmp$chain
  nbchain <- tmp$nbchain
  ordered <- tmp$ordered
  sigsize <- length(siginput)
  rec <- numeric(sigsize)
  plnb <- 0
  
  par(mfrow=c(1,1))
  plot.ts(siginput)
  title("Original signal")
  
  tmp <- matrix(0,nbchain,length(siginput))
  
  totnbnodes <- 0
  idx <- numeric(nbchain)	
  p <- 0
  
  for (j in 1:nbchain){
    phi.x.min <- 2 * 2^(chain[j,3]/nvoice)
    if (chain[j,2] > (para*phi.x.min) ){ 
      
      cat("Chain number",j)
      
      phi.x.max <- 2 * 2^(chain[j,(2+chain[j,2])]/nvoice)
      x.min <- chain[j,1]
      x.max <- chain[j,1] + chain[j,2] - 1
      x.min <- x.min - round(para * phi.x.min)
      x.max <- x.max + round(para * phi.x.max)
      
      tmp2 <- irregrec(siginput[chain[j,1]:(chain[j,1]+chain[j,2]-1)],
                       inputwt[chain[j,1]:(chain[j,1]+chain[j,2]-1),],
                       chain[j,3:(chain[j,2]+2)], compr,noct,nvoice,
                       epsilon, w0 = w0 , prob = prob, fast = fast,
                       para = para, minnbnodes = minnbnodes, real = real);

      totnbnodes <- totnbnodes + tmp2$nbnodes
      
      np <- length(tmp2$sol)
      start <- max(1,x.min)
      end <- min(sigsize,x.min+np-1)
      start1 <- max(1,2-x.min)
      end1 <- min(np, sigsize +1 - x.min)
      end <- end1 - start1 + start
      rec[start:end] <- rec[start:end]+tmp2$sol[start1:end1]
      
      plnb <- plnb + 1
      p <- p + 1
      idx[p] <- j
    }
  }
  
  if(plot == 1) {
    par(mfrow=c(2,1))
    par(cex=1.1)
    plot.ts(siginput)
    title("Original signal")
    title("Reconstructed signal")
  }
  else if (plot == 2) {
    par(mfrow=c(plnb+2,1))
    par(mar=c(2,4,4,4))
    par(cex=1.1)
    par(err=-1)
    plot.ts(siginput)
    title("Original signal")
    
    for (j in 1:p)
      plot.ts(tmp[idx[j],])
    
    plot.ts(Re(rec))	
    title("Reconstructed signal")
   }

  cat("Total number of ridge samples used: ",totnbnodes,"\n")
  par(mfrow=c(1,1))
  list(rec=rec, ordered=ordered, chain = chain, comp=tmp)
}

crcirgrec <- function(siginput, inputgt, beemap, nvoice, freqstep, scale,
                      compr, prob = 0.8, bstep = 5, ptile = 0.01,
                      epsilon = 0.5, fast = TRUE, para = 0,
                      minnbnodes = 3, hflag = FALSE, plot = 2, real = FALSE)
#########################################################################
#     crcirgrec:
#     ----------
#      Reconstruction of a real valued signal from ridges found by 
#      crazy climbers on the gabor transform (irregularly sampled).
#
#      input:
#      ------
#      siginput: input signal
#      inputgt: continuous gabor transform (output of cgt)
#      beemap: output of crazy climber algorithm
#      nvoice: number of different frequencies
#      freqstep: difference between two consecutive frequencies
#      scale: scale of the window
#      compr: subsampling rate for the ridge
#      bstep: used for the chaining
#      ptile: 
#      epsilon: coeff of the Q2 part of reconstruction kernel
#      fast: if set to TRUE, computes the Q2 kernel by Riemann sums
#            if not, uses a Romberg adaptive step quadrature.
#      para:
#      prob: percentile to keep 
#      minnbnodes: minimal number of nodes for a sampled ridge.
#      hflag: if set to FALSE, uses the identity as first term
#             in the reconstruction kernel. If not, uses Q1 instead. 
#      plot: if set to 1, displays the signal, the components, and
#            signal and reconstruction one after another. If set to
#            2, displays the signal, the components, and the
#            reconstruction on the same page. Else, no plot.
#      real: if set to true, uses only real constraints.
#
#      output:
#      -------
#      rec: reconstructed signal
#      ordered: image of the ridges (with different colors)
#      comp: 2D array containing the signals reconstructed from ridges
#
#########################################################################
{
  tmp <- cfamily(beemap,bstep,ptile=ptile)
  image(tmp$ordered)
  chain <- tmp$chain
  nbchain <- tmp$nbchain
  ordered <- tmp$ordered
  sigsize <- length(siginput)
  rec <- numeric(sigsize)
  plnb <- 0
  
  par(mfrow=c(1,1))
  plot.ts(siginput, main="Original signal")
  
  totnbnodes <- 0
  tmp <- matrix(0,nbchain,length(siginput))
  for (j in 1:nbchain){
    if (chain[j,2] > scale){
      nbnodes <- round(chain[j,2]/compr)
      if(nbnodes < minnbnodes)
        nbnodes <- minnbnodes
      
      totnbnodes <- totnbnodes + nbnodes
      
      cat("Chain number", j)
      
      phi.x.min <- scale
      phi.x.max <- scale
      x.min <- chain[j,1]
      x.max <- chain[j,1] + chain[j,2] - 1
      x.min <- x.min - round(para * phi.x.min)
      x.max <- x.max + round(para * phi.x.max)
      x.inc <- 1
      np <- as.integer((x.max-x.min)/x.inc) + 1
      
      tmp2 <- girregrec(siginput[chain[j,1]:(chain[j,1]+chain[j,2]-1)],
                        inputgt[chain[j,1]:(chain[j,1]+chain[j,2]-1),],
                        chain[j,3:(chain[j,2]+2)], nbnodes, nvoice,
                        freqstep, scale, epsilon, fast, prob = prob,
                        para = para, hflag = hflag, real = real)

      start <- max(1,x.min)
      end <- min(sigsize,x.min+np-1)
      start1 <- max(1,2-x.min)
      end1 <- min(np, sigsize +1 - x.min)
      end <- end1 - start1 + start
      rec[start:end] <- rec[start:end]+tmp2$sol[start1:end1]
      tmp[j,start:end] <- tmp2$sol[start1:end1]
      
      plnb <- plnb + 1
      plot.ts(tmp[j,])
    }
    
  }

  if(plot == 1){
    par(mfrow=c(2,1))
    par(cex=1.1)
    plot.ts(siginput, main="Original signal")
    plot.ts(rec, main="Reconstructed signal")
   }
  else if (plot == 2) {
    par(mfrow=c(plnb+2,1))
    par(mar=c(2,4,4,4))
    par(cex=1.1)
    par(err=-1)
    plot.ts(siginput, main="Original signal")
    
    for (j in 1:nbchain)
      plot.ts(tmp[j,])
        
    plot.ts(rec, main="Reconstructed signal")
  }
  
  cat("Total number of ridge samples used: ",totnbnodes,"\n")
  par(mfrow=c(1,1))
  list(rec=rec, ordered=ordered,chain = chain, comp = tmp)
}

