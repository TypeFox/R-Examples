fp.target <-
function(Y=NULL, irfdraws=NULL, nlags=4, constant=TRUE, type="median", labels=colnames(Y), target= TRUE, save=FALSE, legend=TRUE, bands=c(0.16, 0.84), grid=TRUE, bw=FALSE, maxit=1000){ # what else do i need in terms of options?
#
#--- SANITY CHECKS ---#
sanity.check.target(Y=Y, nlags=nlags, irfdraws=irfdraws, constant=constant, type=type, labels=labels, target= target, save=save, legend=legend, bands=bands, grid=grid, bw=bw, maxit=maxit)
#
graphics.off()
par.def <- par(no.readonly = T)
  #--- PARAS---#
  ldg <- legend
  irftype <- type
  goodresp <- irfdraws
  gridgrph <- grid # grid in irf plots 0== none, 1== adds grid to plots
  tgt <- target
  bndtest <- is.null(bands)
  if(bndtest!=TRUE){
    ebupp <- bands[2]# error bands for irf plots
    eblow <- bands[1]
  }else{
    ebupp <- 0.84# error bands for irf plots
    eblow <- 0.16
  }
  varlbl <- labels
  nstep <- dim(irfdraws)[2]
  nvar <- ncol(Y)
  nobs <- nrow(Y)
  n1 <- maxit
  nlags <- nlags
  nnobs0 <- nlags + 1
  nnobs <- nobs - nlags
  nnvar0 <- nvar + 1
#
  if(constant == FALSE){
    CONS <- "F"
    ncoef <- nvar * nlags
    nncoef <- nvar * nlags
    nnvar1 <- nvar * (nlags + 1)
  }else{
    CONS <- "T"
    ncoef <- nvar * (nlags+1)
    nncoef <- nvar * nlags + 1
    nnvar1 <- nvar * (nlags + 1) + 1
  }
 #
  #--- SET UP MATRICES ---#
  targets <- matrix(0.0, nrow=nstep, ncol=nvar)
  rescale <- matrix(0.0,nrow=nstep, ncol=nvar)
  lower  <- matrix(0.0, nrow=nstep, ncol=nvar)
  upper  <- matrix(0.0, nrow=nstep, ncol=nvar)
  #
  #--- GET TARGETS ---#
  for(i in 1:nvar){# that can be done using apply
    for(k in 1:nstep){
      if(irftype=="mean"){
        targets[k,i] <- mean(goodresp[ ,k,i])
      }else{
      targets[k,i] <- quantile(goodresp[ ,k,i], probs=0.5)
      }
      rescale[k,i] <- 1.0/(quantile(goodresp[ ,k,i], probs=0.75) - quantile(goodresp[ ,k,i], probs=0.25))
      lower[k,i]  <- quantile(goodresp[ ,k,i], probs= eblow)
      upper[k,i]  <- quantile(goodresp[ ,k,i], probs= ebupp)
    }
  }
  #
  #--- RFVAR ---#
  model <- rfvar(Y,lags=nlags, const=CONS, breaks=NULL)
  bcoef <- model$By # same order as above but w/const and nvar x nvar x lags
  resid <- model$u # same as above
  data <- model$X
  xx <- model$xx
  #
  #--- SIGMA and SXX ---#
  uu <- crossprod(resid)
  # sigma <- (1/(nnobs-nncoef))*uu
  sigma <- (1/nnobs)*uu
  betaols <- t(bcoef)
  #
  swish <- chol(sigma)
  imfhat <- fn.impulse(betaols, swish, c(nvar, nlags, nstep))
  impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
  #
  #--- PAGAN/FRY GAP ---#
  g <- matrix(1, nrow=nvar-1, ncol=1)
  gapfunc <- minqa::uobyqa(g, fn=FPgap, control = list(maxfun=n1), rescale=rescale, targets=targets, impulses=impulses, nstep= nstep)
  #
  paras <- gapfunc$par
  funcval <- gapfunc$fval
  conv <- gapfunc$ierr
  #
  #--- WARNING MESSAGE ---#
  if(conv>0){
    message('\n Warning! optimisation did not converge.', sep="")
  }
  #
  #--- CREATE MIN IMPULSES ---#
  imp <- matrix(NA, nrow=nstep, ncol=nvar)
  q <- stereo(paras)
  imfhat <- fn.impulse(betaols, swish, c(nvar, nlags, nstep))
  impulses <-  array(imfhat, dim=c(nstep,nvar,nvar))
  for(j in 1:nstep){
    imp[j,] <- t(impulses[j,,]%*%q)
  }
  #
  #--- DETERMINE COLS AND ROWS OF PLOT ---#
  rowsize <-  ceiling(sqrt(nvar))
  colsize <- ceiling(nvar / rowsize)
  #
  #--- GRAPH PARAS ---#
  if(ldg==TRUE){
  par(bty="o", oma = c(4, 1, 1, 1), mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
  }else{
    par(bty="o",  mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
  }
  #--- GENERATE PLOTS ---#
  if(bw==FALSE){
  for(i in 1:nvar){
    ulim <- max(c(imp[,i],upper[,i],lower[,i]))
    llim <- min(c(imp[,i],upper[,i],lower[,i]))
    plot(x=1:nstep, y=imp[,i], type="l", col="red", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
    if(bndtest!=TRUE){
    lines(1:nstep, y=lower[,i], col="blue", lwd=2)
    lines(1:nstep,  y=upper[,i], col="blue", lwd=2)
    }
    abline(h=0, col="black")
    if(target==TRUE){
    lines(1:nstep, y=targets[,i], lty=2, col="blue", lwd=2)
    }
    if(gridgrph==1){
    grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
  }
  }
    if(ldg==TRUE){
      if(target==TRUE){
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if(bndtest==TRUE){
          legend("bottom",  legend=c("Fry-Pagan", "Impulse response"), col = c("red", "blue"), lty=c(1, 2), lwd=c(2, 2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }else{
        legend("bottom",  legend=c("Fry-Pagan", "Impulse response", "Error bands"), col = c("red", "blue", "blue"), lty=c(1, 2, 1), lwd=c(2, 2, 2),horiz=TRUE   , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5 , x.intersp=1, text.width=c(0.25,0.25,0.5))
        }
      }else{
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if(bndtest==TRUE){
        legend("bottom",  legend=c("Fry-Pagan"), col = c("red"), lty=c(1), lwd=c(2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }else{
          legend("bottom",  legend=c("Fry-Pagan",  "Error bands "), col = c("red", "blue"), lty=c(1, 1), lwd=c(2, 2 ), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }
      }
    }
    if(save==TRUE){
      dev.copy(postscript,'fptarget.eps')
      dev.off()
    }
  }else{#start bw
    for(i in 1:nvar){
      ulim <- max(c(imp[,i],upper[,i],lower[,i]))
      llim <- min(c(imp[,i],upper[,i],lower[,i]))
      plot(x=1:nstep, y=imp[,i], type="l", col="black", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
      if(bndtest!=TRUE){
        lines(1:nstep, y=lower[,i],  lty=2, col="black", lwd=2)
        lines(1:nstep,  y=upper[,i],  lty=2, col="black", lwd=2)
      }
      abline(h=0, col="black")
      if(target==TRUE){
        lines(1:nstep, y=targets[,i], col="darkgrey", lwd=2)
      }

      if(gridgrph==1){
        grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
      }
    }
    if(ldg==TRUE){
      if(target==TRUE){
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if(bndtest==TRUE){
          legend("bottom",  legend=c("Fry-Pagan", "Impulse response"), col = c("black", "darkgrey"), lty=c(1, 1), lwd=c(2, 2), horiz=TRUE   , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }else{
          legend("bottom",  legend=c("Fry-Pagan", "Impulse response", "Error bands"), col = c("black", "darkgrey", "black"), lty=c(1, 1, 3), lwd=c(2, 2, 2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5, x.intersp=1, text.width=c(0.25,0.25,0.5))
        }
      }else{
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if(bndtest==TRUE){
          legend("bottom",  legend=c("Fry-Pagan"), col = c("black"), lty=c(1), lwd=c(2), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }else{
          legend("bottom",  legend=c("Fry-Pagan",  "Error bands "), col = c("black", "black"), lty=c(1, 2), lwd=c(2, 2 ), horiz=TRUE , xpd = TRUE,  inset = c(0, 0), bty = "n", cex = 1 , seg.len=4, xjust=0.5)
        }
      }
    }
    if(save==TRUE){
      dev.copy(postscript,'fptarget.eps')
      dev.off()
    }
  }
  par(par.def)
}
