`dldetect` <-
function(input, algorithm, filestem, ...)
{
  chr.names <- names(input$map)
  type <- attr(input, "type")
  n.chr <- length(input$map)
  locations <- list()
  results <- list()
  logfile <- paste(filestem, ".det.log", sep="")
  loc.temp <- list()
  results$converge=TRUE


  if (algorithm=="asreml") {
    n.perm <- input$nperm
    permfile = paste(filestem, ".perm1", sep="")
    write("Note: if n.perm=0, p-values and threshold are analytical; if n.perm>0 these are empirical", logfile)
    write("****************************************************", logfile)
  }

  if (algorithm=="lme") {
    write("Note: p-values and threshold are analytical", logfile)
    write("****************************************************", logfile)
    n.perm <- 0
  }

  test <- dltest(input, algorithm, chrSet=1:n.chr, ...)

  # Fit initial detection step with all chromosomes
  # Output contains p-values for each chromosome variance component
  cat("*******************************************\n")
  cat("* Detection Stage Iteration 1:\n")
  cat("* Testing all chromosomes\n")
  cat("*******************************************\n")

  if (test$converge == FALSE) 	results$converge <- FALSE
  found <- (min(test$adj.pval) < input$alpha)
  chrSet <- which(test$adj.pval < input$alpha)

  if (n.perm>0)
  { 
    write(paste(c("Perm", chr.names), collapse=" "), permfile)
    write(t(cbind(c(1:(n.perm+1)), test$perm.ts)), permfile, append=TRUE)
  }

  # Output observed statistics to logfile
  write(paste("Iteration 1: No. Permutations=",n.perm, sep=""), logfile, append=TRUE)
  write(paste(c("",chr.names), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh, 4), sep=""), logfile, append=TRUE)
  write("Significant chromosomes to be used for scanning/testing:", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)

  itnum = 1
  while (found==1) 
  {
    itnum=itnum+1
    permfile = paste(filestem, ".perm", itnum, sep="")
    loc.temp <- locations

    markers <- vector(length=length(chrSet))
    for (jj in 1:length(chrSet))
    {
	no.qtls <- length(grep(paste("C", chrSet[jj], "M", sep=""), locations))	
	if (type=="f2") no.qtls <- no.qtls/2

  	cat("*******************************************\n")
	cat("Marker Selection Iteration ", itnum, " for ", chr.names[chrSet[jj]],"\n")
  	cat("*******************************************\n")
     	map.results <- dlmapdet(input, algorithm, s.chr=chrSet[jj], chrSet=chrSet, prevLoc=loc.temp, ...)
	
	if (map.results$converge == FALSE) 	results$converge <- FALSE

	sel.mrk <- which.max(map.results$wald)
	markers[jj] <- sel.mrk
	if (type=="f2")
	 locations <- c(locations, paste("C", chrSet[jj], "M", sel.mrk, c("D", "A"), sep="")) else
	 locations <- c(locations, paste("C", chrSet[jj], "M", sel.mrk, sep=""))

	if (itnum==2) {
	  results$profile[[names(input$map)[chrSet[jj]]]] <- rbind(as.vector(input$map[[chrSet[jj]]]), map.results$wald)
	  rownames(results$profile[[names(input$map)[chrSet[jj]]]]) <- c("Position", "Wald")
	}
    } # end of loop over significant chromosomes

    locations <- unlist(locations)

    # Output selected markers to logfile
    write(paste(c("Mrk:", markers), collapse="\t"), logfile, append=TRUE)

    # Do not allow more than one QTL per interval on a chromosome
    ex <- vector()
    chrSet.old <- chrSet

    for (ii in 1:length(chrSet))
 	if ((no.qtls+2)>=length(input$map[[chrSet[ii]]]))
	  ex <- c(ex, ii)

    chrSet <- setdiff(chrSet, chrSet[ex])
  
    if (length(ex)>0)
    {
	write("*******************************************************", logfile, append=TRUE)
	write("Chromosomes full - One QTL located per interval already", logfile, append=TRUE)
	write(paste("Chromosome(s): ", paste(c("", chr.names[chrSet.old[ex]]), collapse="\t"), " removed from testing set", sep=""), logfile, append=TRUE)
	write("*******************************************************", logfile, append=TRUE)
    }

    if (length(chrSet)>0)
    {
	# Test previous subset of chromosomes for significance 
  	cat("*******************************************\n")
	cat("Testing chromosomes for Iteration ", itnum, "\n")
  	cat("*******************************************\n")
  	test <- dltest(input, algorithm, chrSet, prevLoc=locations, ...)

	if (n.perm>0)
	{
  	 write(paste(c("Perm", chr.names[chrSet]), collapse=" "), permfile)
  	 write(t(cbind(c(1:(n.perm+1)), test$perm.ts)), permfile, append=TRUE)
	}

	if (test$converge==FALSE) 	results$converge <- FALSE
  	found <- (min(test$adj.pval) < input$alpha)
  	chrSet <- chrSet[which(test$adj.pval < input$alpha)]
    }

  write("***************************************************************", logfile, append=TRUE)
  write(paste("Iteration ", itnum,": No. Permutations=",n.perm, sep=""), logfile, append=TRUE)

  write("Chromosomes from previous iteration: ", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet.old]), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh,4), sep=""), logfile, append=TRUE)

  if (length(chrSet)>0)
  {
   write("Significant chromosomes for next round of testing/scanning:", logfile, append=TRUE)
   write(paste(c("",chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)
  }
 }  # End of iterative search

  results$qtls$pos <- locations
  return(results)
}

