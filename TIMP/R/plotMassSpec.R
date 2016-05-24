"plotMassSpec" <-
function(multimodel, t, plotoptions)
{
  m <- multimodel@modellist
  resultlist <- multimodel@fit@resultlist 
  
  kinspecerr <- plotoptions@kinspecerr
  superimpose <- plotoptions@superimpose 
  
  if(length(superimpose) < 1 || any(superimpose > length(m)))
    superimpose <- 1:length(m)
  
  allx2 <- allx <- vector() 
  for(i in 1:length(m)) {
    allx2 <- append(allx2, m[[i]]@x2) 
    allx <- append(allx, m[[i]]@x)
  }
  specList <- list() 
  maxs <- mins <- maxspecdim <- 0
  specList <- getSpecList(multimodel, t)
  
  maxA <- lapply(specList, max)
  maxs <- max(unlist(maxA)[superimpose])
  
  minA <- lapply(specList, min)
  mins <- min(unlist(minA)[superimpose])
  
  maxspecdim <- max(unlist(lapply(specList, ncol))[superimpose])

  if(dev.cur() != 1)
    dev.new()
  par(mgp = c(2, 1, 0), mar=c(0,0,0,0), oma = c(1,0,4,0), 
      mfrow=c( ceiling( ncol(specList[[1]]) / 2), 2))		      
  ## make 1 row 2 col if only plotting 2 comp.
  if(ceiling( ncol(specList[[1]]) / 2)==1) par(mfrow=c(2,1))
  xlim <- c(min(allx2),max(allx2))
  ylim <- vector()
  
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
  if(kinspecerr)   
    errtList <- getSpecList(multimodel, t,  getclperr=TRUE) 
  plotted <- FALSE	
  
  for(i in 1) {  ## note that are just plotting dataset's spectra
    if(i %in% superimpose) {		
      if(kinspecerr) {
        if(plotoptions@writeclperr) 
          write.table(errtList[[i]], 
                      file=paste(plotoptions@makeps,
                        "_std_err_clp_", i, ".txt", sep=""), 
                      quote = FALSE, row.names = m[[i]]@x2)
      }
      if (plotoptions@normspec) {
        sp <-  normdat(specList[[i]])
        
      }
      else 
        sp <- specList[[i]]
      for(j in 1:ncol(sp)) {
        if(kinspecerr)
          plotCI(m[[i]]@x2, sp[,j], 
                 uiw=errtList[[i]][,j], 
                 main = "", xlab = plotoptions@xlab,
                 ylab="amplitude", 
                 lty = 1, ##xlim =xlim,ylim=ylim,
                 col = j, sfrac = 0,  type="l", gap = 0,
                 add = !(i == 1), labels = "")
        else {
          if(plotoptions@normspec) {
                      mm<-max( sp[,j] )
                      mi <- min(sp[,j])
                      ylim <- c(mi, mm)
                    }
          
          names.side <- 0 
          if(j <= 2)
            names.side <- 3
                    if(j >= (ncol(sp) - 1)) {
                      names.side <- 1 	
                      par(mar=c(2,0,0,0))
                    }
          if(names.side != 0)	
            barplot3(sp[,j], col = j, border = j, 
                     ylim = ylim, names.arg = m[[i]]@x2, 
                     names.side = names.side, 
                     names.by = plotoptions@axis.by,
                     axes=FALSE)
          else 
            barplot3(sp[,j], col = j, border = j, 
                     ylim = ylim, axes=FALSE)
        }
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
    par(las = 2)
    
  }
  abline(0,0)
  if (dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript		
    
    dev.print(device=pdev, 
              file=paste(plotoptions@makeps, "_massspec.", 
		plotoptions@output, sep=""))
    
  }
}

