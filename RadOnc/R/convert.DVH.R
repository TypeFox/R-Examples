convert.DVH <- function(..., type=NULL, dose=NULL, volume=NULL, dose.units=NULL) {
	type <- match.arg(type, choices=c(NA, "cumulative", "differential"))
	dose <- match.arg(dose, choices=c(NA, "absolute", "relative"))
	volume <- match.arg(volume, choices=c(NA, "relative", "absolute"))
	dose.units <- match.arg(dose.units, choices=c(NA, "cGy", "Gy"))
	arglist <- c(...)
	arglist <- arglist[unlist(lapply(arglist, function(arg) { ((class(arg)[1] %in% c("DVH", "zDVH")) & validObject(arg)) }))]
	N <- length(arglist)
	if (N <=0) {
		return(NULL)
	}
	
	for (i in 1:N) {
		x <- arglist[[i]]
		if (is.empty(x)) { 
			x@dose.type <- dose
			x@volume.type <- volume
			x@type <- type
			x@dose.units <- dose.units
			arglist[[i]] <- x
			next	
		}
		if ((!is.na(dose)) & (dose != x@dose.type)) {
			if (is.na(x@dose.rx)) {
				warning(paste("Cannot convert DVH (", x@structure.name, ") doses because prescription dose is not specified", sep=""), immediate.=TRUE, call.=FALSE)
				arglist[[i]] <- x
				next	
			}
			if (dose == "absolute") {
				x@doses <- x@doses * x@dose.rx / x@rx.isodose
				x@dose.max <- x@dose.max * x@dose.rx / x@rx.isodose
				x@dose.min <- x@dose.min * x@dose.rx / x@rx.isodose
				x@dose.mean <- x@dose.mean * x@dose.rx / x@rx.isodose
				x@dose.median <- x@dose.median * x@dose.rx / x@rx.isodose
				x@dose.mode <- x@dose.mode * x@dose.rx / x@rx.isodose
				x@dose.STD <- x@dose.STD * x@dose.rx / x@rx.isodose
				x@dose.type <- "absolute"
			}
			else {
				x@doses <- x@doses * x@rx.isodose / x@dose.rx
				x@dose.max <- x@dose.max * x@rx.isodose / x@dose.rx
				x@dose.min <- x@dose.min * x@rx.isodose / x@dose.rx
				x@dose.mean <- x@dose.mean * x@rx.isodose / x@dose.rx
				x@dose.median <- x@dose.median * x@rx.isodose / x@dose.rx
				x@dose.mode <- x@dose.mode * x@rx.isodose / x@dose.rx
				x@dose.STD <- x@dose.STD * x@rx.isodose / x@dose.rx
				x@dose.type <- "relative"
			}
		}
		else {
			dose <- x@dose.type
		}
		if ((!is.na(volume)) & (volume != x@volume.type)) {
			if (volume == "absolute") {
				if (x@structure.volume == 0) {
					print(x@structure.volume)
					warning(paste("Cannot convert DVH (", x@structure.name, ") to 'absolute' volume units, because structure has zero volume", sep=""), immediate.=TRUE, call.=FALSE)
				}
				x@volumes <- x@volumes * x@structure.volume / 100
				x@volume.type <- "absolute"
			}
			else {
				x@volumes <- 100 * x@volumes / x@structure.volume
				x@volume.type <- "relative"
			}
		}
		else {
			volume <- x@volume.type
		}
		if ((!is.na(dose.units)) & (dose.units != x@dose.units)) {
			if (dose.units == "cGy") {
				x@dose.rx <- x@dose.rx * 100
				if (x@dose.type == "absolute") {
					x@doses <- x@doses * 100
					x@dose.max <- x@dose.max * 100
					x@dose.min <- x@dose.min * 100
					x@dose.mean <- x@dose.mean * 100
					x@dose.median <- x@dose.median * 100
					x@dose.mode <- x@dose.mode * 100
					x@dose.STD <- x@dose.STD * 100
				}
			}
			else {
				x@dose.rx <- x@dose.rx / 100
				if (x@dose.type == "absolute") {
					x@doses <- x@doses / 100
					x@dose.max <- x@dose.max / 100
					x@dose.min <- x@dose.min / 100
					x@dose.mean <- x@dose.mean / 100
					x@dose.median <- x@dose.median / 100
					x@dose.mode <- x@dose.mode / 100
					x@dose.STD <- x@dose.STD / 100
				}
			}
			x@dose.units <- dose.units
		}
		if ((!is.na(type)) & (type != x@type)) {
			if (type == "cumulative") {
				temp.doses <- x@doses - diff(c(-x@doses[1], x@doses))/2
				x@doses <- c(temp.doses, (2*x@doses - temp.doses)[length(temp.doses)])
				if (volume == "relative") {
					if (class(x) == "DVH") {
						x@volumes <- diffinv(-x@volumes, xi=100)
					}
					else {
						volumes <- diffinv(-x@volumes, xi=matrix(apply(x@volumes, 2, sum), nrow=1))
						class(volumes) <- c("numeric", "matrix")
						colnames(volumes) <- colnames(x@volumes)
						x@volumes <- volumes
					}
				}
				else {
					if (class(x) == "DVH") {
						x@volumes <- diffinv(-x@volumes, xi=x@structure.volume)
					}
					else {
						volumes <- diffinv(-x@volumes, xi=matrix(apply(x@volumes, 2, sum), nrow=1))
						class(volumes) <- c("numeric", "matrix")
						colnames(volumes) <- colnames(x@volumes)
						x@volumes <- volumes
					}
				}
				x@type <- "cumulative"
			}
			else {
				if (class(x) == "DVH") {
					x@volumes <- -diff(x@volumes)
				}
				else {
					volumes <- -apply(x@volumes, 2, diff)
					class(volumes) <- c("numeric", "matrix")
					colnames(volumes) <- colnames(x@volumes)
					x@volumes <- volumes
				}
				x@doses <- x@doses[1:(length(x@doses)-1)] + diff(x@doses)/2
				x@type <- "differential"
			}
		}
		arglist[[i]] <- x
	}
	if (N == 1) { return(arglist[[1]]) }
	return(arglist)
}