setGeneric("LQE",
	function (x, aB, ...) {
		standardGeneric("LQE") 
	}
)

setMethod("LQE", c("ANY", "missing"),
	function (x, aB, ...) {
		stop("argument 'aB' is missing, with no default")
	}
)

setMethod("LQE", c("numeric", "numeric"),
	function (x, aB, fractions=NULL, N=NULL, dose.units=c("cGy", "Gy")) {
		dose.units <- match.arg(dose.units)
		if (is.null(fractions)) {
			if (is.null(N)) {
				warning("arguments 'fractions' and 'N' are missing, with no default")	
				return(NA)	
			}
			else {
				fractions <- N
				dose.units <- paste(dose.units, "N", sep="")
			}
		}
		if (length(fractions) != 2) {
			warning("argument 'fractions' must be of length two")
			return(NA)
		}
		if (any(fractions <= 0)) {
			warning("argument 'fractions' must specify two positive numeric values")
			return(NA)
		}
		if (length(aB) < 1) {
			warning("argument 'aB' is of zero length")
			return(NA)
		}
		if (aB == 0) {
			warning("argument 'aB' must be non-zero")
			return(NA)
		}
		switch(dose.units,
			cGy = {
				return(x * (0.01 * fractions[1] + aB) / (0.01 * fractions[2] + aB))
			},
			Gy = {
				return(x * (fractions[1] + aB) / (fractions[2] + aB))				
			},
			cGyN = {
				return(suppressWarnings(unlist(lapply(x, 
					function (dose) {
						dose <- as.numeric(polyroot(c(-dose * (0.01 * dose / fractions[1] + aB), aB, 0.01 / fractions[2])))
						return(dose[dose >= 0])
					}
 				))))
			},
			GyN = {
				return(suppressWarnings(unlist(lapply(x, 
					function (dose) {
						dose <- as.numeric(polyroot(c(-dose * (dose / fractions[1] + aB), aB, 1 / fractions[2])))
						return(dose[dose >= 0])
					}
 				))))
			}			
		)
	}
)

setMethod("LQE", c("DVH", "numeric"),
	function (x, aB, fractions=NULL, N=NULL, dose.units=c("cGy", "Gy")) {
		dose.units <- match.arg(dose.units)
		if (is.empty(x)) {
			warning("argument 'x' is an empty DVH")
			return(NA)
		}
		x <- convert.DVH(x, dose="absolute", dose.units=dose.units)
		if (is.null(fractions)) {
			if (is.null(N)) {
				warning("arguments 'fractions' and 'N' are both missing, with no default")	
				return(NA)	
			}
			else {
				fractions <- N
				dose.units <- paste(dose.units, "N", sep="")
			}
		}
		if ((length(fractions) != 1) || (fractions <= 0)) {
			warning("argument 'fractions' must specify a single positive numeric value")	
			return(NA)	
		}
		if (length(aB) < 1) {
			warning("argument 'aB' is of zero length")
			return(NA)
		}
		else if (length(aB) > 1) {
			warning(paste("length of 'aB' exceeds length of 'x', will use single value for aB=",aB[1], sep=""))
			aB <- aB[1]
		}
		if (aB == 0) {
			warning("argument 'aB' must be non-zero")
			return(NA)
		}
		switch(dose.units,
			cGy = {
				if ((x$dose.fx == 0) || (x$dose.fx == x$dose.rx / fractions)) {
					x$dose.fx <- x$dose.rx / fractions
					return(x)
				}				
				LQE.calc <- function (doses) {
					return(doses * (0.01 * doses / x$dose.fx + aB) / (0.01 * fractions * doses / x$dose.rx + aB))
				}		
			},
			Gy = {
				if ((x$dose.fx == 0) || (x$dose.fx == x$dose.rx / fractions)) {
					x$dose.fx <- x$dose.rx / fractions
					return(x)
				}				
				LQE.calc <- function (doses) {
					return(doses * (doses / x$dose.fx + aB) / (fractions * doses / x$dose.rx + aB))
				}		
			},
			cGyN = {
				if ((x$dose.fx == 0) || (x$dose.fx == fractions)) {
					x$dose.fx <- fractions
					return(x)
				}				
				LQE.calc <- function (doses) {
					return(suppressWarnings(unlist(lapply(doses, 
						function (dose) {
							dose <- as.numeric(polyroot(c(-dose * (0.01 * dose / x$dose.fx + aB), aB, 0.01 / fractions)))
							return(dose[dose >= 0])
						}
 					))))
				}		
			},
			GyN = {
				if ((x$dose.fx == 0) || (x$dose.fx == fractions)) {
					x$dose.fx <- fractions
					return(x)
				}		
				LQE.calc <- function (doses) {
					return(suppressWarnings(unlist(lapply(doses, 
						function (dose) {
							dose <- as.numeric(polyroot(c(-dose * (dose / x$dose.fx + aB), aB, 1 / fractions)))
							return(dose[dose >= 0])
						}
 					))))
				}		
			}
		)
		x$doses <- LQE.calc(x$doses)
		x$dose.max <- LQE.calc(x$dose.max)
		x$dose.min <- LQE.calc(x$dose.min)
		x$dose.mean <- LQE.calc(x$dose.mean)
		x$dose.median <- LQE.calc(x$dose.median)
		x$dose.mode <- LQE.calc(x$dose.mode)
		x$dose.STD <- LQE.calc(x$dose.STD)
		x$dose.rx <- LQE.calc(x$dose.rx)
		x$dose.fx <- fractions		
		return(x) 
	}
)

setMethod("LQE", c("DVH.list", "numeric"),
	function (x, aB, fractions=NULL, N=NULL, dose.units=NULL) {
		dose.units <- match.arg(dose.units, choices=c("cGy","Gy"), several.ok=TRUE)
		if (length(dose.units) != length(x)) {
			if (length(dose.units) > 1) {
				warning(paste("length of 'x' and 'dose.units' do not match, will use single value for dose.units=",dose.units[1], sep=""))
			}
			dose.units <- rep(dose.units[1], length(x))
		}
		if (length(aB) != length(x)) {
			if (length(aB) > 1) {
				warning(paste("length of 'x' and 'aB' do not match, will use single value for aB=",aB[1], sep=""))
			}
			aB <- rep(aB[1], length(x))
		}
		if (length(fractions) != length(x)) {
			if (length(fractions) > 1) {
				warning(paste("length of 'x' and 'fractions' do not match, will use single value for fractions=",fractions[1], sep=""))
			}
			else if (length(fractions) < 1) {
				if (length(N) != length(x)) {
					if (length(N) > 1) {
						warning(paste("length of 'x' and 'N' do not match, will use single value for N=",N[1], sep=""))
					}
					else if (length(N) < 1) {
						warning("arguments 'fractions' and 'N' are both missing, with no default")	
						return(NA)	
					}
					N <- rep(N[1], length(x))				
				}
				return(new("DVH.list", mapply(LQE, x, aB=aB, N=N, dose.units=dose.units)))			
			}
			fractions <- rep(fractions[1], length(x))
		}
		return(new("DVH.list", mapply(LQE, x, aB=aB, fractions=fractions, dose.units=dose.units)))
	}
)
