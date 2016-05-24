setMethod("plotter", signature(model = "kin"),function(model, 
                                 multimodel, multitheta, plotoptions) {  
  if (!plotoptions@notraces) 
    plotTracesSuper(multimodel, multitheta, plotoptions)
  if(plotoptions@residplot) {
    if(!plotoptions@FLIM) 
      plotResids(multimodel, multitheta, plotoptions) 
  }
  if(plotoptions@residtraces)
    plotTracesResids(multimodel, multitheta, plotoptions)
  if(plotoptions@writefit || plotoptions@writefitivo)
    writeFit(multimodel, multitheta, plotoptions)
  if (length(plotoptions@breakdown)>0) 
    plotKinBreakDown(multimodel, multitheta,plotoptions)  
  if (plotoptions@FLIM) {
    plotFLIM(multimodel, multitheta, plotoptions)
    return()
  }
  if(dev.cur() != 1)
    dev.new()
  par(mgp = c(2, 1, 0), mar=c(3,3,3,2), oma = c(1,0,4,0), 
      mfrow=c(plotoptions@summaryplotrow, plotoptions@summaryplotcol))
  resultlist <- multimodel@fit@resultlist
  m <- multimodel@modellist
  t <- multitheta
  groups <- multimodel@groups
  allx2 <- allx <- vector() 
  for(i in 1:length(m)) {
    allx2 <- append(allx2, m[[i]]@x2) 
    allx <- append(allx, m[[i]]@x)
  }
  xmax <- max(allx)
  xmin <- min(allx)
  x2max <- max(allx2)
  x2min <- min(allx2)
  conmax <- tauList <- muList <- contoplotList <- list()

    if(plotoptions@adddataimage){
            for(i in 1:length(m)) {
		datarev <- t(apply(m[[i]]@psi.df, 2, rev))
		xr <- rev(m[[i]]@x)
		x2lab <- xr[seq(1,length(xr), length=10)]
		x2at <- m[[i]]@x[seq(1,length(xr), length=10)]
		image(x=m[[i]]@x2, y=m[[i]]@x, z=datarev, ylab = plotoptions@xlab, yaxt="n",main = "Data", xlab = plotoptions@ylab, col=tim.colors())
		axis(2, at=x2at, labels=x2lab)
	    }
     }

  if(m[[i]]@anispec$useparperp) 
    calcAniSignal(m, plotoptions)
  f1<-function(x){x[[1]]}   
  f2<-function(x){x[[2]]}
  ## will plot the first concentration from each dataset
  grtoplot <- vector("list", length(m)) 
  for(i in 1:length(m)) {
    cnt <- 1
    notfound <- TRUE 
    while(notfound) {
      for(j in 1:length(groups[[cnt]])) {
        if(groups[[cnt]][[j]][2] == i) {
          grtoplot[[i]]<-list(groups[[cnt]],j)
          notfound<-FALSE
        }
      }
      cnt<-cnt+1
    }	
  }
  for(i in 1:length(m)) {
    group <- grtoplot[[i]][[1]]
    place <-  grtoplot[[i]][[2]]
    dset <- group[[place]][2]
    irfmu <- unlist(lapply(resultlist[[i]]@irfvec, f1))
    irftau <- unlist(lapply(resultlist[[i]]@irfvec, f2))
    muList[[i]] <- irfmu 
    tauList[[i]] <- irftau
    if(plotoptions@writerawcon || 
       length(plotoptions@writeplaincon)>0){ 
      C <- getKinConcen(group, multimodel, t, doConstr = FALSE, oneDS = place)
      if(plotoptions@writerawcon)
        write.table(C, file=paste(plotoptions@makeps,
                         "_rawconcen_dataset_", dset, ".txt", sep=""),
                    quote = FALSE,
                    row.names = m[[dset]]@x)
      if(length(plotoptions@writeplaincon)>0){
        xplot <- plotoptions@writeplaincon$x 
        CWRITE <- C 
        write.table(CWRITE, file=paste(plotoptions@makeps,
                              "_plaincon_dataset_", dset, ".txt", sep=""),
                    quote = FALSE,
                    row.names = m[[dset]]@x)
      }
    }
    contoplotList[[length(contoplotList)+1]] <- getConToPlot( getKinConcen(
        group, multimodel, t, oneDS = place), m[[i]]@cohspec, m[[i]]@cohcol)
    conmax[[i]] <- attributes(contoplotList[[length(contoplotList)]])$max
  }
  minc <- do.call("min", contoplotList)
  maxc <- do.call("max", contoplotList)
  for(i in 1:length(m)) {
    matlinlogplot(m[[i]]@x, contoplotList[[i]], muList[[i]][1], 
                  plotoptions@linrange, 
                  type = "l", add = !(i == 1), lty = i, ylab = "concentration", 
                  xlab = plotoptions@xlab, main = "Concentrations", 
                  xlim = c(xmin, xmax), ylim = c(minc, maxc) )
    if(plotoptions@writecon) 
      write.table(contoplotList[[i]], file=paste(plotoptions@makeps,
                                        "_concen_dataset_", i, ".txt",
                                        sep=""), quote = FALSE,
                  row.names = m[[i]]@x)
  }
  if(length(model@kin2scal)!=0) {
    perA<-vector()
    for(i in 1:length(m)) 
      perA <- append(perA, t[[i]]@kin2scal[2])
    matplot(1:length(m), perA, type = "l", 
            main = "% Concentration Photoconverted", 
            xlab = "Dataset number", ylab = "percent")
  }
  spectralist <- getSpecList(multimodel, t)
  specList <- list() 
  for (i in 1:length(m)) {
    if (length(conmax) > 0) 
      spec <- getSpecToPlot(spectralist[[i]], conmax[[i]], 
                            m[[i]]@cohcol,  plotoptions@plotcohcolspec)
    else spec <- spectralist[[i]]
    
    specList[[length(specList)+1]] <- spec
  }
  maxs <- do.call("max", specList)
  mins <- do.call("min", specList)
  for (i in 1:length(m)) {
    matplot(m[[i]]@x2, specList[[i]], type = "l", main = "Spectra", 
            xlab = plotoptions@ylab, ylab = "amplitude", lty = i, 
            add = !(i == 1), xlim=c(x2min,x2max), ylim=c(mins, maxs))
    if(plotoptions@writespec)
      write.table(specList[[i]], file=paste(plotoptions@makeps,
                                   "_spec_dataset_", i, ".txt", sep=""),
                  row.names = m[[i]]@x2, quote=FALSE) 
  }
  abline(0,0)
  for (i in 1:length(m)) {
    matplot(m[[i]]@x2, normdat(specList[[i]]), type = "l", 
	    main = "Normalized spectra",  xlim=c(x2min,x2max),
            xlab = plotoptions@ylab, ylab = "amplitude", lty = i, 
            add = !(i == 1))
    abline(0,0)
    if(plotoptions@writenormspec)
      write.table( normdat(specList[[i]]),	
                  file=paste(plotoptions@makeps,
                    "_normspec_dataset_", i, ".txt", sep=""),
                  row.names = m[[i]]@x2, quote=FALSE)  
  }
  abline(0,0)
  notplotted <- TRUE
  minmu <- do.call("min", muList)
  maxmu <- do.call("max", muList)
  for (i in 1:length(m)) {
    if(m[[i]]@dispmu) {
      matplot(m[[i]]@x2, muList[[i]], type = "l", 
              main = "IRF location",  xlim=c(x2min,x2max),
              add = !notplotted, ylim = c(minmu, maxmu), col=i,
              xlab = plotoptions@ylab, ylab = "IRF location")
      notplotted <- FALSE  
    }
  }
  notplotted <- TRUE
  mintau <- do.call("min", tauList)
  maxtau <- do.call("max", tauList)
  for (i in 1:length(m)) {
    if(m[[i]]@disptau) {
      
      matplot(m[[i]]@x2, tauList[[i]], type = "l", 
              main = "IRF width",   add = !notplotted, col=i,
              xlim=c(x2min,x2max),ylim = c(mintau, maxtau),
              xlab = plotoptions@ylab, ylab = "IRF width")
            notplotted <- FALSE
    }
  }
  svddatalist <- list()
  residlist <- svdresidlist <- list()
  for (i in 1:length(m)) {
    residuals <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
    for (j in 1:length(resultlist[[i]]@resid)) {
      residuals[, j] <- resultlist[[i]]@resid[[j]]
    }
    svdresidlist[[i]] <- doSVD(residuals, 2, 2)
    residlist[[i]] <- residuals
    svddatalist[[i]] <- doSVD(multimodel@data[[i]]@psi.df, 2, 2)
  }
  
  rangeleftd <- rangerightd <- rangevald <- 0
  rangeleftr <- rangerightr <- rangevalr <- 0
  for(i in 1:length(m)) {
    rangeleftr <- range(svdresidlist[[i]]$right[1,], rangeleftr)
    rangerightr <- range(svdresidlist[[i]]$right[1,], rangerightr)
    rangevalr <- range(svdresidlist[[i]]$values, rangevalr)

    rangeleftd <- range(svddatalist[[i]]$right[1,], rangeleftd)
    rangerightd <- range(svddatalist[[i]]$right[1,], rangerightd)
    rangevald <- range(svddatalist[[i]]$values, rangevald)
  }
  ##START RESID PLOTTING
  for (i in 1:length(m)) {          
    limd<- max(  max(residlist[[i]]), abs(min(residlist[[i]]))) 
    if (! (any(diff(m[[i]]@x) <= 0) || any(diff(m[[i]]@x2) <= 0)))
      image.plot(m[[i]]@x, m[[i]]@x2, 
                 residlist[[i]], xlab = plotoptions@xlab, 
                 ylab = plotoptions@ylab, 
                 main = paste("Residuals Dataset", i),  
                 zlim=c(-limd,limd),
                 col = diverge_hcl(40, h = c(0, 120), c = 60, 
                   l = c(45, 90), power = 1.2))
  }
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
      xpos <- m[[i]]@x 
      xpos[which(xpos <= 0)] <- NA  
      matplot(xpos, svdresidlist[[i]]$left[,1], 
              type = "l", ylim = rangeleftr, 
              main = "LSV Residuals",  add = !(i == 1),
              log = "x", xlab = plotoptions@xlab, col=i, 
              ylab = plotoptions@ylab)
    }
  }
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
      matplot(m[[i]]@x2, svdresidlist[[i]]$right[1,], 
              type = "l", xlim=c(x2min,x2max),
              ylim = rangerightr, 
              main = "RSV Residuals", 
              xlab = plotoptions@xlab, add = !(i == 1),
              col=i, ylab = plotoptions@ylab)
    }
  }
  for (i in 1:length(m)) {     
    if(i == 1)
      plot(1:length(svdresidlist[[i]]$values), 
           log10(svdresidlist[[1]]$values), xlab="",
           ylab = plotoptions@ylab, col=i, ylim = rangevalr,
           main = "SVal Residuals", type = "b")
    else
      lines(1:length(svdresidlist[[i]]$values), 
            log10(svdresidlist[[i]]$values), type = "b", col=i)
  }
  ##START DATA PLOTTING
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
      xpos <- m[[i]]@x 
      xpos[which(xpos <= 0)] <- NA  
      matplot(xpos, svddatalist[[i]]$left[,1], type = "l", 
              main = "LSV Data",  add = !(i == 1),
              log = "x", xlab = plotoptions@xlab, col=i, 
              ylim = rangeleftd, 
              ylab = plotoptions@ylab)
    }
  }    
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
	    matplot(m[[i]]@x2, svddatalist[[i]]$right[1,], type = "l", 
                    main = "RSV Data", xlim=c(x2min,x2max),
                    ylim = rangerightd, 
                    xlab = plotoptions@xlab, add = !(i == 1),
                    col=i, ylab = plotoptions@ylab)
          }
  }
  for (i in 1:length(m)) {     
    if(i == 1)
      plot(1:length(svddatalist[[i]]$values), 
           log10(svddatalist[[1]]$values), xlab="",
           ylab = plotoptions@ylab, col=i, ylim=rangevald, 
           main = "SingVal Data", type = "b")
    else
      lines(1:length(svddatalist[[i]]$values), 
            log10(svddatalist[[i]]$values), type = "b", col=i)
  }
  ## OTHER PLOTS
  for (i in 1:length(m)) {
    pl <- FALSE
    if (length(m[[i]]@dscalspec$perclp) != 0) 
      if (m[[i]]@dscalspec$perclp) {
        if (!pl) {
          plot(m[[i]]@x2, t[[i]]@drel, 
               main = "Dataset scaling per clp", 
               xlab = plotoptions@ylab, ylab = "", type = "l")
          pl <- TRUE
        }
        else lines(m[[i]]@x2, t[[i]]@drel, type = "l", 
                   col = i)
      }
  }
  if(length(plotoptions@title) != 0){
    tit <- plotoptions@title
    if(plotoptions@addfilename) tit <- paste(tit,m[[i]]@datafile)
  }
  else {
    tit <- ""
    if(plotoptions@addfilename) tit <- paste(tit, m[[i]]@datafile)
  }
  mtext(tit, side = 3, outer = TRUE, line = 1)
  
  writeEst(multimodel, multitheta, plotoptions)
  # TODO: replace functionality provided by displayEst
  displayEst(plotoptions)
  
  if (dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf
    else  pdev <- postscript
    dev.print(device = pdev,
              file = paste(plotoptions@makeps, 
                "_summary.",  plotoptions@output,
                sep = ""))
  }
  if (plotoptions@plotkinspec) {
    plotClp(multimodel, t, plotoptions)
  }
  if (plotoptions@kinspecest) {
    plotKinSpecEst(t, plotoptions, multimodel)
  }
  
}
          )
