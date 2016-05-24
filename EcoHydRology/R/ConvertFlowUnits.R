ConvertFlowUnits <-
function(cfs=NULL, cmd=NULL, cms=NULL, WA, mmd=NULL, AREAunits = "mi2"){
	feetpermile <- 5280
	mmperfoot <- 304.8
	secondsperday <- 3600 * 24
	kmpermile <- 1.609  # 
	if (AREAunits == "mi2" &  (!is.null(cmd) | !is.null(cms))) {
		WA <- WA * kmpermile^2
		AREAunits <- "km2"
	}
	if (!is.null(cfs) & AREAunits == "mi2"){
		mm_d <- cfs * secondsperday * mmperfoot / WA / (feetpermile^2) 
		return(mm_d)
	} else if (! is.null(cmd)){  
		mm_d <- cmd /WA/1000
		return(mm_d)
	} else if (! is.null(cms)){ 
		mm_d <- cms * secondsperday/WA/1000
		return(mm_d)
	} else { 
		cfs <- mmd * WA * (feetpermile^2) / (secondsperday * mmperfoot)
		return(cfs)
	}
}
