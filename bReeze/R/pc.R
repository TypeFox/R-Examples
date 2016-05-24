pc <- function(pc, ...) {
	# this is not really the way the S3 system should work, is it?
	if(is.list(pc) || is.data.frame(pc)) r <- pc.default(pc, ...)
	else if(is.character(pc)) r <- pc.read(pc, ...)
	else stop(substitute(pc) , "must be a numeric vector of wind speeds, the name of, or the path to a 'wgt' or 'pow' file containing power curve data")
	return(r)
}


pc.default <-
function(pc, rho=1.225, rated.p, desc, ...) {
###	creating power curve object

	v <- p <- cp <- ct <- NULL
	if(any(names(pc)!=c(v, p, cp, ct))) stop("Variable(s) not recognized - power curve data may contain 'v', 'p', 'cp' and/or 'ct'")
	if(length(names(pc))>0) {
		if(any(names(pc)=="")) stop("Please name all power curve variables")
		if(any(names(pc)=="v")) v <- pc$v
		else stop("Wind speed 'v' is mandatory")
		if(any(names(pc)=="p")) p <- pc$p
		else stop("Power 'p' is mandatory")
		if(any(names(pc)=="cp")) cp <- pc$cp
		if(any(names(pc)=="ct")) ct <- pc$ct
	} else {
		if(is.list(pc)) {
			if(length(pc)>0) v <- pc[[1]]
			else stop("Wind speed 'v' is mandatory")
			if(length(pc)>1) p <- pc[[2]]
			else stop("Power 'p' is mandatory")
			if(length(pc)>2) cp <- pc[[3]]
			if(length(pc)>3) ct <- pc[[4]]
		} else if(is.data.frame(pc)) {
			if(length(pc)>0) v <- pc[,1]
			else stop("Wind speed 'v' is mandatory")
			if(length(pc)>1) p <- pc[,2]
			else stop("Power 'p' is mandatory")
			if(length(pc)>2) cp <- pc[,3]
			if(length(pc)>3) ct <- pc[,4]
		}
	}
	if(!is.vector(v)) stop("'v' requires numeric vector")
	if(!is.vector(p)) stop("'p' requires numeric vector")
	if(length(v)!=length(p)) stop("Different vector length of 'v' and 'p'")
	if(!is.null(cp)) if(!is.vector(cp)) stop("'cp' requires numeric vector")
	if(!is.null(cp)) if(length(v)!=length(cp)) stop("Different vector length of 'v' and 'cp'")
	if(!is.null(ct)) if(!is.vector(ct)) stop("'ct' requires numeric vector")
	if(!is.null(ct)) if(length(v)!=length(ct)) stop("Different vector length of 'v' and 'ct'")
	if(missing(rho)) rho <- 1.225
	if(missing(rated.p)) rated.p <- max(p, na.rm=TRUE)
	if(missing(desc)) desc <- NULL
	
	pc <- data.frame(cbind(v, p, cp, ct))
	names <- c("v", "P")
	attr(pc, "units") <- c("m/s", "kW")
	if(!is.null(cp)) {
		names <- append(names, "cp")
		attr(pc, "units") <- append(attr(pc, "units"), "-")
	}
	if(!is.null(ct)) {
		names <- append(names, "ct")
		attr(pc, "units") <- append(attr(pc, "units"), "-")
	}
	names(pc) <- names
	attr(pc, "rho") <- rho
	attr(pc, "rated.power") <- rated.p
	if(!is.null(desc)) attr(pc, "description") <- desc
	attr(pc, "call") <- list(func="pc.default")
	
	class(pc) <- "pc"
	return(pc)
}


pc.read <-
function(pc, ...) {
### importing power curve from WAsP .wgt file or WindPower program .pow file
		
	if(system.file(package="bReeze", "powercurves", pc)!="") pc <- system.file(package="bReeze", "powercurves", pc)
	if(!file.exists(pc)) stop("File not found")
	
	type <- substr(pc, nchar(pc)-3, nchar(pc))
	if(!any(c(".pow", ".wtg")==type)) stop("Cannot handle file - only WAsP .wtg files and WindPower program .pow files are supported")
	
	r <- NULL
	if(type==".pow") {
		pow <- read.table(pc, as.is=TRUE)
		cut.out <- as.numeric(pow[4,1])
		if(is.na(cut.out) || is.null(cut.out)) stop("Cannot handle file")
		v <- seq(1, cut.out, 1)
		p <- tail(pow, -1)
		suppressWarnings(if(is.na(as.numeric(tail(p, 1)))) p <- head(p, -1))
		p <- as.numeric(p[5:(cut.out+4),1])
		desc <- pow[1,1]
		r <- pc.default(pc=list(v=v, p=p), rho=1.225, desc=desc)
		attr(r, "call") <- list(func="pc.read", pc=pc)
	} else if(type==".wtg") {
		stopifnot(requireNamespace("XML", quietly=TRUE))
		wtg <- XML::xmlTreeParse(pc, asTree=TRUE)
		if(is.null(wtg$doc$children$WindTurbineGenerator)) stop("Cannot handle file")
		n <- length(wtg$doc$children$WindTurbineGenerator)
		idx <- 3
		if(n>4) {
			rho <- NULL
			for(i in 3:(n-1)) rho <- append(rho, as.numeric(XML::xmlAttrs(wtg$doc$children$WindTurbineGenerator[[i]])[["AirDensity"]]))
			idx <- which.min(abs(rho-1.225))+2
		}
		rho <- as.numeric(XML::xmlAttrs(wtg$doc$children$WindTurbineGenerator[[idx]])[["AirDensity"]])
		n <- length(wtg$doc$children$WindTurbineGenerator[[idx]][["DataTable"]])
		v <- p <- ct <- NULL
		for(i in 1:n) {
			v <- append(v, as.numeric(XML::xmlAttrs(wtg$doc$children$WindTurbineGenerator[[idx]][["DataTable"]][[i]])[["WindSpeed"]]))
			p <- append(p, as.numeric(XML::xmlAttrs(wtg$doc$children$WindTurbineGenerator[[idx]][["DataTable"]][[i]])[["PowerOutput"]])/1000)
			ct <- append(ct, as.numeric(XML::xmlAttrs(wtg$doc$children$WindTurbineGenerator[[idx]][["DataTable"]][[i]])[["ThrustCoEfficient"]]))
		}
		desc <- XML::xmlAttrs(XML::xmlRoot(wtg))[["Description"]]
		r <- pc.default(pc=list(v=v, p=p, ct=ct), rho=rho, desc=desc)
		attr(r, "call") <- list(func="pc.read", pc=pc)
	}
	
	class(r) <- "pc"	
	return(r)
}
