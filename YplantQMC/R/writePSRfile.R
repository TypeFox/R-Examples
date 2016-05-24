
writePSRfile <- function(yplantsim, suffix="", ndigits=6){

	if(!inherits(yplantsim, "yplantsim"))
		stop("Need object output by YplantDay.")
	
	psr <- yplantsim$psrdata
	pfile <- yplantsim$plant$pfile
	lfile <- yplantsim$plant$lfile
	
	proot <- gsub("\\.p$","",pfile, ignore.case=TRUE)
  if(is.null(proot))
    proot <- format(Sys.time(), "%Y-%m%-%d")
  
	filen <- paste0(proot,suffix,"-YplantQMC.PSR")
	
	r <- c()
	r[1] <-  "Simulation summary for plant constructed from files:"
	r[2] <-  paste("   P file:", pfile)
	r[3] <-  paste("   L file:", lfile)
	r[4] <- "----------------------------------------------------\n"
	totAPAR <- sum(psr$PARleaf * psr$LAplant * psr$timestep) * 10^-6
	totdifAPAR <- sum(psr$PARdiff * psr$LAplant * psr$timestep) * 10^-6
	totdirAPAR <- sum(psr$PARdir * psr$LAplant * psr$timestep) * 10^-6
	
	if("A" %in% names(psr)){
	totA <- sum(psr$A * psr$LAplant * psr$timestep) * 10^-6
	r[5] <- paste("Total photosynthesis      (mol day-1)  :", round(totA,3))
	} else r[5] <- ""
	r[6] <- paste("Total absorbed PAR (APAR) (mol day-1)  :", round(totAPAR,3))
	r[7] <- paste("Total Diffuse APAR        (mol day-1)  :", round(totdifAPAR,3))
	r[8] <- paste("Total Direct APAR         (mol day-1)  :", round(totdirAPAR,3))
	r[9] <- "----------------------------------------------------\n"
	r[10] <- "Variables by timestep:"
	r[11] <- "(Some may not be included)"
	r[12] <- "timeofday   - Time of day (hours) of timestep."
	r[13] <- "PAR0        - Incident PAR above canopy (mu mol m-2 s-1)"
	r[14] <- "PARinc      - Incident PAR below canopy (mu mol m-2 s-1)"
	r[15] <- "PARleaf     - Absorbed total PAR (mu mol m-2 s-1)"
	r[16] <- "PARdir      - Absorbed direct PAR (mu mol m-2 s-1)"
	r[17] <- "PARdiff     - Absorbed diffuse PAR (mu mol m-2 s-1)"
	r[18] <- "A           - CO2 assimilation rate (mu mol m-2 s-1)"
    r[19] <- "A0          - A by horizontal leaf above canopy (mu mol m-2 s-1)"
    r[20] <- "E           - Transpiration rate (mmol m-2 s-1)"  
    r[21] <- "timestep    - Length of timestep (seconds)" 
	r[22] <- "LAplant     - Total plant leaf area (m2)"   
    r[23] <- "LAproj      - Total projected leaf area (m2) ('NA' if diffuse PAR only)."
    r[24] <- "LAsunlit    - Total displayed (sunlit) leaf area (m2) ('NA' if diffuse PAR only)."
	r[25] <- ""
	
	writeLines(r, filen)
	options(warn=-1)
	write.table(round(psr,ndigits), filen, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE)
	options(warn=0)
	message("Results written to file: ", filen)

return(invisible(filen))
}


