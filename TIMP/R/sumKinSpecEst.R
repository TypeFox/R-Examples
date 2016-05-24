"sumKinSpecEst"  <- function (listFits, addtitle=TRUE, 
customtitle="", preps="", ylimlist=list(), kinspecerr=TRUE ) {

  if(dev.cur() != 1)
    dev.new()
  par(mfrow=c(length(listFits),2))
  par(mai = c(0.5, 0.6, .5, 0.5))
  par(mgp = c(2, 1, 0), mar=c(1,3,1,1))
  par(oma = c(1,0,4,0))
  
  if(length(customtitle) != 0) 
    tit <- customtitle 
  else
    tit <- "Summary of EADS for: "
	max_x2 <- max(listFits[[1]]$currModel@modellist[[1]]@x2)
  min_x2 <- min(listFits[[1]]$currModel@modellist[[1]]@x2)
  limlist <- list()
  for(i in 1:length(listFits))
    limlist[[i]] <- vector()			  
  if(length(ylimlist) > 0) 
    for(i in 1:length(ylimlist))
      limlist[[ylimlist[[i]]$ind]] <- ylimlist[[i]]$ylim
  for(i in 1:length(listFits)){
    m <- listFits[[i]]$currModel@modellist
    for(j in 1:length(m)){
      max_x2 <- max(m[[j]]@x2, max_x2)
      min_x2 <- min(m[[j]]@x2, min_x2)
    }
  }
  for(i in 1:length(listFits)){
    resultlist <- listFits[[i]]$currModel@fit@resultlist 
    multimodel <- listFits[[i]]$currModel
    m <- multimodel@modellist
    t <- listFits[[i]]$currTheta
	      plotoptions <- listFits[[i]]$currModel@optlist[[1]]
    if(length(plotoptions@addest) == 0)
      plotoptions@addest <- c("kinpar")
    
    plotClp(multimodel, t, plotoptions, kinspecerr=kinspecerr,
            newplot=FALSE, min_clp=min_x2,max_clp=max_x2, ylim=limlist[[i]])

    ##plotEst(multimodel, plotoptions) 
    
    if(length(customtitle)==0) {
      if(addtitle && i==1)
		  tit <- paste(tit, "  ", plotoptions@title, sep="")
      else if(addtitle) 
        tit <- paste(tit, plotoptions@title, sep=", ")
    }
  }
  if (addtitle) {
    mtext(tit, side = 3, outer = TRUE, line = 1)
    par(las = 2)
  }
	if (dev.interactive() && length(plotoptions@makeps) != 0) {
          dev.print(device = postscript, file = paste(preps, 
                                           "_kinspecsum.ps", sep = ""))
        }
}
