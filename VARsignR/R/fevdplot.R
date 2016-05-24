fevdplot <-
function(fevddraws=NULL, type="median", labels=unlist(dimnames(fevddraws)[3]), save=FALSE, bands=c(0.16, 0.84), grid=TRUE, bw=FALSE, table=FALSE, periods=NULL){

    #--- SANITY CHECK ---#
 # sanity.check.fevdplot(fevddraws=fevddraws, type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw, table=table)


    graphics.off()
  par.def <- par(no.readonly = T)
    #graph parameter
    goodresp <- fevddraws
    irftype <- type #  0== median, 1== mean response
    gridgrph <- grid # grid in irf plots 0== none, 1== adds grid to plots
    bndtest <- is.null(bands)
    if(bndtest!=TRUE){
    ebupp <- bands[2]# error bands for irf plots
    eblow <- bands[1]
    }else{
      ebupp <- 0.84# error bands for irf plots
      eblow <- 0.16
    }
    varlbl <- labels
    nstep <- dim(fevddraws)[2]
    nvar <- dim(fevddraws)[3]
    periodst <- is.null(periods)

    if(irftype=="mean"){
      imp_responses <- array(NA, dim=c(3, nstep, nvar))
      irfbands <- apply(goodresp,c(2,3),quantile,probs=c(eblow, ebupp))
      irfmean <-  array(apply(goodresp,c(2,3),mean), dim=c(1,nstep, nvar))
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
      imp_responses[1,,] <- irfmean
      imp_responses[2:3,,] <- irfbands
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
    }else{
      imp_responses <- apply(goodresp,c(2,3),quantile,probs=c(0.5, eblow, ebupp))
      dimnames(imp_responses) <- list(c("FEVD", "lower", "upper"),1:nstep, varlbl)
    }
    impt <- imp_responses
    impt <- aperm(impt,c(3,2,1))

    if(table==FALSE){
    #--- DETERMINE COLS AND ROWS OF PLOT ---#
    rowsize <-  ceiling(sqrt(nvar))
    colsize <- ceiling(nvar / rowsize)

    #-- GENERATE PLOTS ---#
    # dev.off()
    par(bty="o", mfcol=c(rowsize, colsize), mar=c(rep(2.5,4)))
    if(bw==FALSE){
      for(i in 1:nvar){
        ulim <- max(impt[i,,1:3])
        llim <- min(impt[i,,1:3])
        plot(x=1:nstep, y=impt[i,,1], type="l", col="red", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
        lines(1:nstep, impt[i,,2], col="blue", lwd=2)
        lines(1:nstep, impt[i,,3], col="blue", lwd=2)
        }
        abline(h=0, col="black")
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
      }
      if(save==TRUE){
        dev.copy(postscript,'fevd.eps')
        dev.off()
      }
    }else{
      for(i in 1:nvar){
        ulim <- max(impt[i,,1:3])
        llim <- min(impt[i,,1:3])
        plot(x=1:nstep, y=impt[i,,1], type="l", col="black", lwd=2, xlab="", ylab="", main= paste(varlbl[i]), ylim=c(llim, ulim))
        if(bndtest!=TRUE){
        lines(1:nstep, impt[i,,2], lty=2, col="black", lwd=2)
        lines(1:nstep, impt[i,,3], lty=2, col="black", lwd=2)
        }
        abline(h=0, col="black")
        if(gridgrph==1){
          grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",   lwd = par("lwd"))
        }
      }
      if(save==TRUE){
        dev.copy(postscript,'fevd.eps')
        dev.off()
      }

    }
    par(par.def)
    }else{#if table ==TRUE
    if(periodst==1){# if period ==null
      fevdtable <- t(impt[,,1])
      fevdtable <- round(fevdtable, 2)
      colnames(fevdtable) <- varlbl
    }else{# if period !=null
      fevdtable <- t(impt[,periods,1])
      fevdtable <- round(fevdtable, 2)
      colnames(fevdtable) <- varlbl
    }
      #returns table
  return(fevdtable)
  }
}
