"plotSpectraSuper" <- function(multimodel, multitheta, plotoptions)
{
  if(dev.cur() != 1)
    dev.new()
  m <- multimodel@modellist   
  t <- multitheta   
  res <- multimodel@fit@resultlist
  superimpose <- plotoptions@superimpose 
  if(length(superimpose) < 1 || any(superimpose > length(m)))
    superimpose <- 1:length(m)
  divdrel <- plotoptions@divdrel
  allx2 <- allx <- vector() 
  for(i in superimpose) {
    allx2 <- append(allx2, m[[i]]@x2) 
    allx <- append(allx, m[[i]]@x)
  }
  if(length(plotoptions@xlim) == 0) {
    x2max <- max(allx2)
    x2min <- min(allx2)
  }
  allx <- sort(unique(allx))
  ## set up plot layout 
  par(plotoptions@paropt)
  par(mgp = c(2, 1, 0), mar=c(0,2,3,0), oma = c(1,0,4,0))
  if(length(plotoptions@selectedtraces) > 0 ) {
    seltraces <- plotoptions@selectedtraces 
    xx <- vector()
    for(i in superimpose) {
      xx <- append(m[[i]]@x[seltraces],xx)
    }
    lensel <- length(unique(xx))
  }	
  else {
    seltraces <- 1:length(allx) 
    lensel <- length(seltraces)
  }
  ymax <- ymin <- rep(0,length(allx)) # is including 0 good?
  
  # ugly, fix

  for (j in 1:length(allx)) {
    for(i in 1:length(m)) {
      k <- which(m[[i]]@x  == allx[j])[1]
      if(i %in% superimpose && k %in% seltraces) {
        data <-  m[[i]]@psi.df[k,]
        fitted <- res[[i]]@fitted[[k]]
        if(divdrel && length(t[[i]]@drel)!=0){
          if(length(m[[i]]@dscalspec$perclp)!=0) 
            if(m[[i]]@dscalspec$perclp)	{	
              data <- data/t[[i]]@drel[k]
              fitted <- fitted/t[[i]]@drel[k]
            }
            else {
              data <- data/t[[i]]@drel
              fitted <- fitted/t[[i]]@drel
            }
        }
        if(m[[i]]@weight) 
          fitted <- fitted/m[[i]]@weightM[, k]
        ymax[j] <- max(fitted,data,ymax[j])
        ymin[j] <- min(fitted,data,ymin[j])
        
      }
    }
  }
  par(mfrow = n2mfrow(lensel))
  for (j in 1:length(allx)) {
    plotted <- FALSE 
    for(i in 1:length(m)) {
      k <- which(m[[i]]@x == allx[j])[1]
      
      if(i %in% superimpose && k %in% seltraces) {	
        data <- m[[i]]@psi.df[k,] 
        fitted <- res[[i]]@fitted[[k]] 
        if(divdrel && length(t[[i]]@drel)!=0) 
          if(length(m[[i]]@dscalspec$perclp)!=0) 
            if(m[[i]]@dscalspec$perclp)	{	
              data <- data/t[[i]]@drel[k]
              fitted <- fitted/t[[i]]@drel[k]
            }
            else {
              data <- data/t[[i]]@drel
              fitted <- fitted/t[[i]]@drel
            }
        if(m[[i]]@weight) 
          fitted <- fitted/m[[i]]@weightM[k,]
        cold <- if(length(plotoptions@coldata) == 0) i
         else plotoptions@coldata[i]
        ltd <- if(length(plotoptions@ltydata) == 0) 
          1 else plotoptions@ltydata 
        colf <- if(length(plotoptions@colfit) == 0)
          if(length(m) > 1) i else 2 else plotoptions@colfit[i]
        ltf <- if(length(plotoptions@ltyfit) == 0) 
          3 else plotoptions@ltyfit
        if(!plotted){
          plot(m[[i]]@x2, data, 
               type = "l", 
               ylim=c(ymin[j],ymax[j]), 
               xlab = plotoptions@xlab,  
               ylab ="amplitude", col = cold, lty = ltd)
          title(signif(m[[i]]@x[k]), line=0)
          plotted <- TRUE
        }
        else 
          lines(m[[i]]@x2, data, type = "l", col = cold, lty = ltd)
        lines(m[[i]]@x2, fitted, type = "l", col = colf, lty = ltf)
      }
    }
  }
  if(length(plotoptions@title) != 0){
    tit <- plotoptions@title
    if(plotoptions@addfilename)
      tit <- paste(tit,m[[i]]@datafile)
  }
  else {
    tit <- ""
    if(plotoptions@addfilename) tit <- paste(tit, m[[i]]@datafile)
  }
  mtext(tit, side = 3, outer = TRUE, line = 1)
  par(las = 2)
  
                                        # MAKE PS
  if(dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript
    dev.print(device=pdev, 
              file=paste(plotoptions@makeps, "_selectedtraces.", 
		plotoptions@output,
		sep=""))
  }
}

