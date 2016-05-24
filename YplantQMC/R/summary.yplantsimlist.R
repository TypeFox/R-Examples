#'@method summary yplantsimlist
#'@S3method summary yplantsimlist
summary.yplantsimlist <- function(object, totals=TRUE, writefile=FALSE, ...){

	x <- object
  
	# Sum these variables.
	vars <- c("PARleaf","PARdir","PARdiff","PARinc","A","A0","E")
	
  # Toss the ones that are not in the psrdata (e.g, E is not always present)
  vars <- intersect(vars, names(x[[1]]$psrdata))
	
  if(totals){
  	sums <- lapply(x, function(z){
    		psr <- z$psrdata
    		v <- colSums(psr[,vars] * psr$timestep * psr$LAplant) *10^-6
    		names(v) <- paste0("tot",vars)
    		return(v)
  	})
    
  	sumdfr <- as.data.frame(do.call("rbind",sums))	
  	pfdfr <- data.frame(pfile=sapply(x, function(x)x$plant$pfile),
  	                    lfile=sapply(x, function(x)x$plant$lfile))
    sumdfr <- cbind(pfdfr, sumdfr)
    
  	if(writefile){
  	
  		filen <- paste0("YplantDayBatchResults-",as.Date(Sys.time()),".txt")
  		unlink(filen)
  	
  		r <- c()
  		r[1] <- "pfile - plant file used to construct plant."
  		r[2] <- "lfile - leaf file used to construct plant."
  		r[3] <- "totPARleaf - total PAR absorption (mol day-1)" 
  		r[4] <- "totPARdfir - total direct PAR absorption (mol day-1)"
  		r[5] <- "totPARdiff - total diffuse PAR absorption (mol day-1)"
      r[6] <- "totPARinc  - total incident PAR (mol day-1)"
      r[7] <- "totA - total CO2 assimilation (mol day-1)"
  		r[8] <- "totA0 - total CO2 assimilation by unshaded horizontal leaf (mol day-1)"
  		r[9] <- "totE - total transpiration (mmol day-1)"
  		
  		r <- c("YplantDay Bath Simulation Result - produced with YplantQMC\n\n",r)
  		r <- c(r, "\n\n")
  	}
  }
  
	if(!totals){
	  sums <- lapply(x, function(z){
	    psr <- z$psrdata
      v <- colMeans(psr[,vars])
	    names(v) <- paste0("mean",vars)
	    return(v)
	  })
	  
	  sumdfr <- as.data.frame(do.call("rbind",sums))	
	  pfdfr <- data.frame(pfile=sapply(x, function(x)x$plant$pfile),
	                      lfile=sapply(x, function(x)x$plant$lfile))
	  sumdfr <- cbind(pfdfr, sumdfr)
	  
	  if(writefile){
	    
	    filen <- paste0("YplantDayBatchResults-",as.Date(Sys.time()),".txt")
	    unlink(filen)
	    
	    r <- c()
	    r[1] <- "pfile - plant file used to construct plant."
	    r[2] <- "lfile - leaf file used to construct plant."
	    r[3] <- "meanPARleaf - mean PAR absorption (mumol m-2 s-1)" 
	    r[4] <- "meanPARdfir - mean direct PAR absorption (mumol m-2 s-1)"
	    r[5] <- "meanPARdiff - mean diffuse PAR absorption (mumol m-2 s-1)"
      r[6] <- "meanPARinc  - mean incident PAR (mumol m-2 s-1)"
	    r[7] <- "meanA - mean CO2 assimilation (mumol m-2 s-1)"
	    r[8] <- "meanA0 - mean CO2 assimilation by unshaded horizontal leaf (mumol m-2 s-1)"
	    r[9] <- "meanE - mean transpiration (mmol m-2 s-1)"
	    
	    r <- c("YplantDay Bath Simulation Result - produced with YplantQMC\n\n",r)
	    r <- c(r, "\n\n")
	  }
	}

  
  if(writefile){
		writeLines(r, filen)
		
		options(warn=-1)
		write.table(sumdfr, filen, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE)
		options(warn=0)

		message("\nSimulation results (plant totals) written to file:\n ", filen)
		return(invisible(sumdfr))
	} else {
		return(sumdfr)
	}
}