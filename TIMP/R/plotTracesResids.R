"plotTracesResids" <- function(multimodel, multitheta, plotoptions)
{
  if(dev.cur() != 1)
    dev.new()
  m <- multimodel@modellist   
  t <- multitheta   
  res <- multimodel@fit@resultlist
  superimpose <- plotoptions@superimpose 
  if(length(superimpose) < 1 || any(superimpose > length(m)))
    superimpose <- 1:length(m)  
  allx2 <- allx <- vector() 
  for(i in superimpose) {
    allx2 <- append(allx2, m[[i]]@x2) 
    allx <- append(allx, m[[i]]@x)
  }
  if(length(plotoptions@xlim) == 0) {
    xmax <- max(allx)
    xmin <- min(allx)
  } 
  allx2 <- sort(unique(allx2))
  ## set up plot layout 
  par(plotoptions@paropt)
  par(mgp = c(2, 1, 0), mar=c(0,2,3,0), oma = c(1,0,4,0))
  if(length(plotoptions@selectedtraces) > 0 ) {
    seltraces <- plotoptions@selectedtraces 
    xx <- vector()
	for(i in superimpose) {
          xx <- append(m[[i]]@x2[seltraces],xx)
	}
	lensel <- length(unique(xx))
  }	
  else {
    seltraces <- 1:length(allx2) 
    lensel <- length(seltraces)
  }
  resList <- vector("list",length=length(m))
  maxr <- vector(length=length(m))
  for(i in 1:length(m)){
    resList[[i]] <- as.matrix(unlist(res[[i]]@resid),
                                     byrow=(m[[i]]@clpType=="x"))
    dim(resList[[i]]) <- dim(m[[i]]@psi.df) 
  }
  xx <- do.call(rbind,resList)
  ymin <- apply(xx, 2, min)
  ymax <- apply(xx, 2, max)
  par(mfrow = n2mfrow(lensel))
  for (j in 1:length(allx2)) {
    plotted <- FALSE 
    for(i in 1:length(m)) {
      k <- which(m[[i]]@x2 == allx2[j])[1]
      if(i %in% superimpose && k %in% seltraces) {	        
        irfmu <- res[[i]]@irfvec[[k]][1]		
        rr <- resList[[i]][,k]
        if(!plotted){
          linlogplot(m[[i]]@x, rr, irfmu, 
                     plotoptions@linrange, type="l", 
                     ylim=c(ymin[j],ymax[j]),xlim=c(xmin,xmax), 
                     xlab = plotoptions@xlab,  lty=1,
                     ylab ="amplitude", col = i)
          title(signif(m[[i]]@x2[k]), line=0)
          plotted <- TRUE
        }
        else
          lines(linloglines(m[[i]]@x, irfmu, plotoptions@linrange), 
		       rr, col = i, type = "l", lty = 1)
      }
    }
  }
  if(length(plotoptions@title) != 0){
    tit <- plotoptions@title
    if(plotoptions@addfilename) tit <- paste(tit,m[[i]]@datafile)
  }
  else {
    tit <- ""
    if(plotoptions@addfilename)
      tit <- paste(tit, m[[i]]@datafile)
  }
  mtext(tit, side = 3, outer = TRUE, line = 1)
  par(las = 2)
  ## MAKE PS
  if(dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript
    dev.print(device=pdev, 
              file=paste(plotoptions@makeps, "_tracesresids.", 
		plotoptions@output,
		sep=""))
  }
}

