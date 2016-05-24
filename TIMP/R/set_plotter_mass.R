setMethod("plotter", signature(model = "mass"), function(model, 
                                 multimodel, multitheta, plotoptions) {
  if(length(plotoptions@paropt) == 0)
    plotoptions@paropt <- par(mgp = c(1.5, 1, 0), 
			      mai = c(0.5, 0.6, .5, 0.5),mar=c(3,3,4,1))
  plotoptions@addest <- c("peakpar")
  resultlist <- multimodel@fit@resultlist
  if(plotoptions@nummaxtraces > 0 ) 
    plotoptions@selectedtraces <- getSelectedTracesMax(multimodel, 
                                                       multitheta, plotoptions)
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
    plotKinBreakDown(multimodel, multitheta, plotoptions)           
  if(dev.cur() != 1)
    dev.new()
  par(mgp = c(2, 1, 0), mar=c(3,3,3,2), oma = c(1,0,4,0), 
      mfrow=c(plotoptions@summaryplotrow, plotoptions@summaryplotcol))
  m <- multimodel@modellist
  t <- multitheta
  allx2 <- allx <- vector() 
  for(i in 1:length(m)) {
    allx2 <- append(allx2, m[[i]]@x2) 
    allx <- append(allx, m[[i]]@x)
  }
  xmax <- max(allx)
  xmin <- min(allx)
  x2max <- max(allx2)
  x2min <- min(allx2)
  
  conList <- list() 
  minc <- maxc <- 0 
  for (i in 1:length(m)) {
    C <- compModelMass(theta = t[[i]], model = m[[i]])
    if(plotoptions@writerawcon) 
      write.table(C, file=paste(plotoptions@makeps,
                       "_rawconcen_dataset_", i, ".txt", sep=""),
                  quote = FALSE,
                  row.names = m[[i]]@x)
    CTil <- C
    
    conList[[length(conList)+1]] <- CTil
    minc <- min(minc, min(conList[[length(conList)]]))
    maxc <- max(maxc, max(conList[[length(conList)]]))
  }
  if(plotoptions@scale.concen) {
    specList <- getSpecList(multimodel, t)	
    conList <- scaleConList(conList, specList) 
  }
  for(i in 1:length(m)) {
    matlinlogplot(m[[i]]@x, conList[[i]], 0, 
                  plotoptions@linrange, 
                  type = "l", add = !(i == 1), lty = i, ylab = "concentration", 
                  xlab = plotoptions@xlab, main = "Concentrations", 
                  xlim = c(xmin, xmax), ylim = c(minc, maxc) )
    if(plotoptions@writecon) 
      write.table(conList[[i]], file=paste(plotoptions@makeps,
             "_concen_dataset_", i, ".txt", sep=""), quote = FALSE,
                  row.names = m[[i]]@x)
  }
  svddatalist <- list()
  residlist <- svdresidlist <- list()
  for (i in 1:length(m)) {
    residuals <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
    for (j in 1:length(resultlist[[i]]@resid)) {
      residuals[, j] <- resultlist[[i]]@resid[[j]]
    }
    svdresidlist[[length(svdresidlist) + 1]] <- doSVD(residuals, 2, 2)
    residlist[[length(residlist) + 1]] <- residuals
    svddatalist[[length(svddatalist) + 1]] <- doSVD(multimodel@data[[i]]@psi.df, 2, 2)
        }
  maxleftr <- minleftr <- maxrightr <- minrightr <- maxvalr <- minvalr <- 0
  maxleftd <- minleftd <- maxrightd <- minrightd <-  maxvald <- minvald <- 0
  for(i in 1:length(m)) {
    maxleftr <- max(svdresidlist[[i]]$left[,1], maxleftr)
	      maxrightr <- max(svdresidlist[[i]]$right[1,], maxrightr)
    minleftr <- min(svdresidlist[[i]]$left[,1], minleftr)
    minrightr <- min(svdresidlist[[i]]$right[1,], minrightr)
    maxvalr <- max(svdresidlist[[i]]$values, maxvalr)
    minvalr <- min(svdresidlist[[i]]$values, minvalr)
    
    maxleftd <- max(svddatalist[[i]]$left[,1], maxleftd)
    maxrightd <- max(svddatalist[[i]]$right[1,], maxrightd)
    minleftd <- min(svddatalist[[i]]$left[,1], minleftd)
    minrightd <- min(svddatalist[[i]]$right[1,], minrightd)
    maxvald <- max(svddatalist[[i]]$values, maxvald)
    minvald <- min(svddatalist[[i]]$values, minvald)
    
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
              type = "l", ylim = c(minleftr, maxleftr), 
              main = "Left sing. vectors residuals ",  add = !(i == 1),
              log = "x", xlab = plotoptions@xlab, col=i, 
              ylab = plotoptions@ylab)
      
    }
  }
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
      
      matplot(m[[i]]@x2, svdresidlist[[i]]$right[1,],
              type = "l", xlim=c(x2min,x2max),
              ylim = c(minrightr, maxrightr), 
              main = "Right sing. vectors residuals ", 
              xlab = plotoptions@xlab, add = !(i == 1),
              col=i, ylab = plotoptions@ylab)
    }
  }
  for (i in 1:length(m)) {     
	    if(i == 1)
	      plot(1:length(svdresidlist[[i]]$values), 
                   log10(svdresidlist[[1]]$values), xlab="",
                   ylab = plotoptions@ylab, col=i, 
                   main = "Sing. values residuals", type = "b")
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
                    main = "Left sing. vectors data",  add = !(i == 1),
                    log = "x", xlab = plotoptions@xlab, col=i, 
                    ylim = c(minleftd, maxleftd),
                    ylab = plotoptions@ylab)
          }
  }    
  for (i in 1:length(m)) {     
    if (m[[i]]@nt > 1 && m[[i]]@nl > 1) {
      matplot(m[[i]]@x2, svddatalist[[i]]$right[1,], type = "l", 
              main = "Right sing. vectors data", xlim=c(x2min,x2max),
              ylim = c(minrightd, maxrightd),
              xlab = plotoptions@xlab, add = !(i == 1),
              col=i, ylab = plotoptions@ylab)
    }
  }
       for (i in 1:length(m)) {     
         if(i == 1)
           plot(1:length(svddatalist[[i]]$values), 
	      log10(svddatalist[[1]]$values), xlab="",
                ylab = plotoptions@ylab, col=i,
                main = "Sing. values data", type = "b")
         else
           lines(1:length(svddatalist[[i]]$values), 
	      log10(svddatalist[[i]]$values), type = "b", col=i)
       }
  if (dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript
    dev.print(device = pdev, file = paste(plotoptions@makeps, 
                                       "_summary.",  
                               plotoptions@output,
                               sep = ""))
  }       	
  plotMassSpec(multimodel, t, plotoptions)
})
