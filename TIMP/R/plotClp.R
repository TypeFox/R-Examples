"plotClp" <-
function(multimodel, t, plotoptions, newplot=TRUE, max_clp=NA, min_clp=NA, 
ylim=vector(), kinspecerr=FALSE)
{
  m <- multimodel@modellist
  resultlist <- multimodel@fit@resultlist 
  if(newplot) {
    par(oma = c(3,3,3,0), mfrow= c(1,1))
    kinspecerr <- plotoptions@kinspecerr
    if(dev.cur() != 1)
      dev.new()
  }
  superimpose <- plotoptions@superimpose 
  if(length(superimpose) < 1 || any(superimpose > length(m)))
            superimpose <- 1:length(m)
  clp <- m[[1]]@clpType
  if(clp == "x") inp <- "x2" else inp <- "x" 
  if(is.na(max_clp) || is.na(max_clp))
    withlim <- FALSE 
  else		 withlim <- TRUE
  allclp <- allinp <- vector() 
  for(i in 1:length(m)) {
    allclp <- append(allclp, slot(m[[i]], clp) )
    allinp <- append(allinp, slot(m[[i]], inp))
  }
  specList <- getSpecList(multimodel, t)
  
  if(kinspecerr)   
    errtList <- getSpecList(multimodel, t,  getclperr=TRUE) 
  for(i in 1:length(m)) {
    cohcol <- m[[i]]@cohcol 
    spec <- getSpecToPlot(specList[[i]], 1, 
                          cohcol)
    
    if(!identical(m[[i]]@cohcol, 0)) {	      
      spec <- spec[,-cohcol]
      
      if(kinspecerr)
        errtList[[i]] <-  errtList[[i]][,-cohcol]
    }
    specList[[i]] <- spec
	
  }	      			      
  maxA <- lapply(specList, max)
  maxs <- max(unlist(maxA)[superimpose])
  
  minA <- lapply(specList, min)
  mins <- min(unlist(minA)[superimpose])

  maxspecdim <- max(unlist(lapply(specList, ncol))[superimpose])
  if(!withlim) 
    xlim <- c(min(allclp),max(allclp))
  else xlim <- c(min_clp, max_clp)
  if(length(plotoptions@xlimspec) == 2) 
    xlim <- plotoptions@xlimspec
  if(length(plotoptions@ylimspec) == 2) 
    ylim <- plotoptions@ylimspec
  if(length(ylim) == 0) 
    ylim <- c(mins, maxs)
  if(length(plotoptions@ylimspecplus) == 2)
    ylim <- ylim + plotoptions@ylimspecplus
  if (plotoptions@normspec)  
    ylim <- c(-1,1)
  
  plotted <- FALSE
  
  for(i in 1:length(m)) {
    if(i %in% superimpose) {		
      if(kinspecerr) {
        if(plotoptions@writeclperr) 
          write.table(errtList[[i]], 
                      file=paste(plotoptions@makeps,
                        "_std_err_clp_", i, ".txt", sep=""), 
                      quote = FALSE, row.names = slot(m[[i]], clp))
        
      }
      if (plotoptions@normspec) 
        sp <-  normdat(specList[[i]])
      else 
        sp <- specList[[i]]
      for(j in 1:ncol(sp)) {
        if(plotoptions@specinterpol) { 
          xx <- predict(interpSpline(slot(m[[i]],clp), sp[,j],
                                     bSpline=plotoptions@specinterpolbspline),
			nseg = plotoptions@specinterpolseg)
          
          if(!plotted) {
            plot(xx, lty = if(plotoptions@samespecline) 1
            else i, main = "", xlab = plotoptions@xlab,
                 ylab="", xlim =xlim,ylim=ylim, col = j,
                 pch= if(plotoptions@specinterpolpoints) 26-j else NA,
                 type=if(plotoptions@specinterpolpoints) "b" else
                 "l")
          }
          else lines(xx, col = j, 
                     lty = if(plotoptions@samespecline) 1 else i,
                     pch= if(plotoptions@specinterpolpoints) 26-j else NA,
                     type=if(plotoptions@specinterpolpoints) "b" else
                     "l")
          if(kinspecerr)
            plotCI(slot(m[[i]],clp),
                   sp[,j], uiw=errtList[[i]][,j], pch
                   = if(plotoptions@specinterpolpoints) 26-j else NA,
                   col = j, sfrac = 0, type="p", gap = 0, add =TRUE,
                   labels = "", lty = if(plotoptions@samespecline) 1
                   else i)
          if(plotoptions@writespecinterpol) 
            write.table(xx$y, file=paste(plotoptions@makeps,
                                "dataset", i, "_smoothedspectracomponent_",
                                j, ".txt", sep=""), 
                        quote = FALSE, row.names = xx$x)
        }
        else 
          if(kinspecerr)
            plotCI(slot(m[[i]],clp), sp[,j], 
                   uiw=errtList[[i]][,j], 
                   main = "", xlab = plotoptions@xlab,
                   ylab="", lty = if(plotoptions@samespecline) 1
                   else i, xlim =xlim,ylim=ylim,
                   col = j, sfrac = 0,  type="l", gap = 0,
                   add = !(i == 1 &&	j == 1), labels = "")
          else {
            if(!plotted) {
              plot(slot(m[[i]],clp), sp[, j], 
                   lty = if(plotoptions@samespecline) 1 else i, 
                   main = "", xlab = plotoptions@xlab,
                   ylab="", xlim =xlim,ylim=ylim, col = j,
                   type="l")
              plotted<-TRUE 
            }
            else lines(slot(m[[i]],clp), sp[, j], col = j, 
                       lty = if(plotoptions@samespecline) 1 else i)
            
          }
        plotted <- TRUE 
      }
    }
  }
  if (newplot && length(plotoptions@title) != 0) {
	  if(length(plotoptions@title) != 0){
            tit <- plotoptions@title
            if(plotoptions@addfilename) tit <- paste(tit,m[[i]]@datafile)
          }
          else {
            tit <- ""
            if(plotoptions@addfilename) tit <- paste(tit, m[[i]]@datafile)
          }
          mtext(tit, side = 3, outer = TRUE, line = 1)
          par(las = 2)
    	  
        }
  abline(0,0)
  if (dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript		
    
    dev.print(device=pdev, 
              file=paste(plotoptions@makeps, "_kinspec.", 
		plotoptions@output, sep=""))

  }
}

