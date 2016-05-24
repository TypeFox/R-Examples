#'@method summary plant3dlist
#'@S3method summary plant3dlist
#'@rdname summary.plant3d
summary.plant3dlist <- function(object, writefile=FALSE, ...){

	plants <- object
	nplants <- attributes(plants)$nplants
	sums <- list()
	message("Summarizing ",nplants," plants.\n")
	message("Overall progress:")
	wp <- txtProgressBar(title = "Summarizing plants.", label = "", 
        min = 0, max = nplants, initial = 0, width = 50, style=3)
	for(i in 1:nplants){
		sums[[i]] <- try(summary(plants[[i]], quiet=TRUE, progressbar=FALSE, ...))
		if(inherits(sums[[i]], "try-error")){close(wp);stop("Problem with plant number ",i)}
		setTxtProgressBar(wp, i)
	}
	close(wp)

	sums2 <- lapply(1:length(sums[[1]]), function(x)sapply(sums, "[[", x))
	sumsdfr <- as.data.frame(sums2)
	names(sumsdfr) <- names(sums[[1]])
	sumsdfr$nsignif <- NULL

	
	if(!writefile)return(sumsdfr)
	
	if(writefile){
	
		filen <- paste0("PlantSummaries-",as.Date(Sys.time()),".txt")
		unlink(filen)
		
		sumsdfr$ALAC <- with(sumsdfr, LA / crownsurf)
		
		sumsdfr$stemsurf <- sumsdfr$stemsurf * 10^4
		sumsdfr$stemvol <- sumsdfr$stemvol * 10^6
	
		r <- c()
		r[1] <- "  LA           -  Total leaf area (m2)"
		r[2] <- "  meanleafsize -  Mean leaf size (cm2)"
		r[3] <- "  nleavesp     -  Number of leaves"
		r[4] <- "  leaflen      -  Mean leaf length (cm)"
		r[5] <- "  meanleafang  -  Mean leaf angle (deg)"
		r[6] <- "  wmeanleafang -  Mean leaf angle weighted by leaf area (deg)"
		r[7] <- "  Xellipsoid   -  Ellipsoidal leaf angle dist. par."
		r[8] <- "  crownvol     -  Crown volume (convex hull) (m3)"
		r[9] <- "  crownsurf    -  Crown surface area (convex hull) (m2)"
		r[10] <- "  ALAC         -  Crown density (AL/AC) (m2 m-2)"
		r[11] <- "  cw           -  Crown width (m)"
		r[12] <- "  cl           -  Crown length (m)"
		r[13] <- "  htot         -  Total height (m)"
		r[14] <- "  cshape       -  Crown shape index (-)"
		r[15] <- "  stemsurf     -  Stem + branch surface area (cm2)"
		r[16] <- "  stemvol      -  Stem + branch volume (cm3)"
		r[17] <- "  stemdiam     -  Stem base diameter (mm)"
		r[18] <- "  meanpath     -  Mean pipe length (mm)"
		r[19] <- "  sdpath       -  Standard deviation of pipe length (mm)"
		r[20] <- "  totlen       -  Total woody segment length (mm)"
		r[21] <- "  Ek           -  Expected distance to 5 nearest leaves (no edge corr.)"
		r[22] <- "  Ek2          -  Expected distance to 5 nearest leaves (with edge corr.)"
		r[23] <- "  Ok           -  Observed distance to 5 nearest leaves"
		r[24] <- "  disp         -  Dispersion parameter (no edge corr.)"
		r[25] <- "  disp2        -  Dispersion parameter (with edge corr.)"
		if(!is.na(sums[[1]]$STARbar))
		     r[26] <- "  STARbar      -  Spherically averaged STAR (m2 m-2)"
		
		r <- c("Plant summary results  - produced with YplantQMC.\n", r)
		r <- c(r, "\n\n\n")
		writeLines(r, filen)
		
		options(warn=-1)
		write.table(sumsdfr, filen, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE)
		options(warn=0)

		message("\nSummary written to file:\n ", filen)
		
		return(invisible(sumsdfr))
	}

}