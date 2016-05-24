

calcparhrly <- function(zens, DOY, totpar, KHRS=24, TAU=0.76, 
	fbeam_daily=NA, fbeam_method=c("spitters","constant")){

	fbeam_method <- match.arg(fbeam_method)
	
	# Rarely, missing values (solar disk visible, center of solar disk below horizon?)
	zens[is.na(zens)] <- pi/2
	
	# Beam fraction for daily PAR (Spitters 1986 method).
	if(is.na(fbeam_daily) & fbeam_method == "spitters"){
		f <- .Fortran("CALCFBMD",
                  as.integer(DOY),
                  as.double(zens),
                  as.double(totpar),
                  as.double(-999),
                  as.integer(KHRS),
                  PACKAGE="YplantQMC")
		fbeam_daily <- f[[4]]
	}

	# fbeam varies throughout the day; following Spitters et al (1986) algorithm
	
	# Simulated PAR data.
	partot_beam <- fbeam_daily * totpar
	partot_diff <- (1-fbeam_daily)* totpar

	fbeams <- rep(-999, KHRS)
	PARs <- rep(-999, KHRS)

	f2 <- .Fortran("CALCPARHRLY",
                 as.double(partot_beam),
                 as.double(partot_diff),
                 as.double(zens),
                 as.double(PARs),
                 as.double(fbeams),
                 as.double(TAU),
                 as.integer(KHRS),
                 PACKAGE="YplantQMC")

	PAR0 <- f2[[4]]
	fbeam <- f2[[5]]

	# fbeam is constant throughout the day (input).
	if(fbeam_method == "constant"){
		fbeam <- rep(fbeam_daily, length(zens))
	}
	return(data.frame(PAR0=PAR0, fbeam=fbeam))
}


