"plotSpecKin" <-
function(multimodel, t, plotoptions, newplot=TRUE, max_x=NA, min_x=NA, 
ylim=vector(), kinspecerr=FALSE)
{
	m <- multimodel@modellist
	resultlist <- multimodel@fit@resultlist 
	if(newplot) {
	   par(oma = c(0,0,2,0))
	   kinspecerr <- plotoptions@kinspecerr
	}
	superimpose <- plotoptions@superimpose 
        if(length(superimpose) < 1 || any(superimpose > length(m)))
           superimpose <- 1:length(m)
	if(is.na(max_x) || is.na(max_x))
			 withlim <- FALSE 
	else		 withlim <- TRUE
        allx2 <- allx <- vector() 
	for(i in 1:length(m)) {
	  allx2 <- append(allx2, m[[i]]@x2) 
	  allx <- append(allx, m[[i]]@x)
	}
	specList <- list() 
	maxs <- mins <- maxspecdim <- 0
	specList <- getSpecList(multimodel, t)
	for(i in 1:length(m)) {
	   spec <- specList[[i]]
	   if(i %in% superimpose) { 	
		      maxs <- max(maxs, max(spec))
		      mins <- min(mins, min(spec))
		      maxspecdim <- max(maxspecdim, ncol(spec))
	   }
        }	      			      
	if(!withlim) 
		xlim <- c(min(allx),max(allx))
	else xlim <- c(min_x, max_x)
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
        for(i in 1:length(m)) {
             if(i %in% superimpose) {		
		if(kinspecerr) {
			if(plotoptions@writeclperr) 
		          write.table(errtList[[i]], 
			  file=paste(plotoptions@makeps,
		          "_std_err_clp_", i, ".txt", sep=""), 
		          quote = FALSE, row.names = m[[i]]@x)
	        }
		if (plotoptions@normspec) 
				    sp <-  normdat(specList[[i]])
		      else 
				    sp <- specList[[i]]
	      
	      for(j in 1:ncol(sp)) {
		    if(plotoptions@specinterpol) { 
		       xx <- predict(interpSpline(m[[i]]@x, 
			sp[,j], bSpline=plotoptions@specinterpolbspline),
			nseg = plotoptions@specinterpolseg)
			
		       if(!plotted) {
			plot(xx, lty = if(plotoptions@samespecline) 1
			else i, main = "", xlab = plotoptions@xlab,
			ylab="amplitude", xlim =xlim,ylim=ylim, col = j,
			type="l")

			}
			else lines(xx, col = j, 
			lty = if(plotoptions@samespecline) 1 else i)
			if(kinspecerr)
			 plotCI(m[[i]]@x, sp[,j], uiw=errtList[[i]][,j], pch
			 = if(plotoptions@specinterpolpoints) 26-j else NA,
			 col = j, sfrac = 0, type="p", gap = 0, add =TRUE,
			 labels = "", lty = if(plotoptions@samespecline) 1
			else i)
		     
		      if(plotoptions@writespecinterpol) 
		          write.table(xx$y, file=paste(plotoptions@makeps,
		          "dataset", i, 
			  "_smoothedspectracomponent_", j, ".txt", sep=""), 
		          quote = FALSE, row.names = xx$x)
		     }
		     else 
		      if(kinspecerr)
		       plotCI(m[[i]]@x, sp[,j], 
		       uiw=errtList[[i]][,j], 
		       main = "", xlab = plotoptions@xlab,
		       ylab="amplitude", lty = if(plotoptions@samespecline) 1
			else i, xlim =xlim,ylim=ylim,
		       col = j, sfrac = 0,  type="l", gap = 0,
		       add = !(i == 1 &&	j == 1), labels = "")
		      else {
		       if(!plotted) {
			plot(m[[i]]@x, sp[, j], 
			lty = if(plotoptions@samespecline) 1 else i, 
			main = "", xlab = plotoptions@xlab,
			ylab="amplitude", xlim =xlim,ylim=ylim, col = j,
			type="l")
                        plotted<-TRUE 
		       }
			else lines(m[[i]]@x, sp[, j], col = j, 
			lty = if(plotoptions@samespecline) 1 else i)

		     }
		     plotted<-TRUE 
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
            dev.print(device = pdev, file = paste(plotoptions@makeps, 
                "_speckin.", plotoptions@output, sep = ""))
        }
}

