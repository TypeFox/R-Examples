`dllocalize` <-
function(input, algorithm, QTLperChr, ...)
{
  locations <- list()
  results <- list()
  type <- attr(input, "type")
  dfMerged <- input$dfMerged
  
  chrSet <- which(QTLperChr>0)

  # Initialize convergence flag
  results$converge=TRUE

  # Loop over all chromsomes which were significant in the preceding detection step
    for (kk in 1:length(chrSet))
    {
     	no.qtls <- QTLperChr[chrSet[kk]][[1]]
	
	mrk <- grep(paste("C", chrSet[kk], "M", sep=""), names(dfMerged))
	chr <- sort(c(mrk, grep(paste("C", chrSet[kk], "P", sep=""), names(dfMerged))))

	if (type=="f2") {
	  mrk <- mrk[seq(1, length(mrk), 2)]
	  chr <- chr[seq(1, length(chr), 2)]
	}

  	cat("*******************************************\n")
 	cat("Scanning chromosome ", names(input$map)[chrSet[kk]], " for ", no.qtls, " QTL\n") 
 	cat("*******************************************\n")
	
 	# note: loc1 only stores info for each chromosome, locations stores it for the entire thing
	loc1 <- list()
	# 1D scan
	for (mm in 1:no.qtls)
	{
	  # 1D scan for QTL
     	  map.results <- dlmaploc(input, algorithm, s.chr=chrSet[kk], chrSet=chrSet, prevLoc=loc1, ...)
	
	  # Update convergence flag
	  if (map.results$converge == FALSE) 	results$converge <- FALSE

	  # Find position of maximum wald score and marker markers
	  # or alternately use the likelihood for localization

	  sel.pos <- which.max(map.results$wald)
	#  if (type=="f2") sel.pos <- sel.pos*2-1

	  # mark.l/r are now the index of the marker which is on the left or right of the selected position
	  if (sel.pos==1) mark.l <- mrk[1] else 
		mark.l <- max(mrk[mrk<chr[sel.pos]])
	  if (sel.pos==length(chr)) mark.r <- mrk[length(mrk)] else
	  	mark.r <- min(mrk[mrk>chr[sel.pos]])

	  if (type=="f2")
	    loc1$mrk <- c(loc1$mrk, names(dfMerged)[c(mark.l:(mark.l+1), mark.r:(mark.r+1))]) else
	  loc1$mrk <- c(loc1$mrk, names(dfMerged)[c(mark.l, mark.r)])

	  if (type=="f2")
	    loc1$pos <- c(loc1$pos, names(dfMerged)[chr[sel.pos]:(chr[sel.pos]+1)]) else
	  loc1$pos <- c(loc1$pos, names(dfMerged)[chr[sel.pos]])

	  if (mm==1)
	  {
		results$profile[[names(input$map)[chrSet[kk]]]] <- t(rbind(input$mapp[[chrSet[kk]]], map.results$wald))
		colnames(results$profile[[names(input$map)[chrSet[kk]]]]) <- c("Position", "Wald")
	  }
  	}

	locations$mrk <- c(locations$mrk, loc1$mrk)
	locations$pos <- c(locations$pos, loc1$pos)

    } # end of loop over significant chromosomes
 
  results$qtls <- locations
  return(results)
}

