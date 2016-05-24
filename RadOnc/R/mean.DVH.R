setGeneric("mean",
	mean
)

setMethod("mean", "DVH.list",
	function (x, ..., type=c("cumulative", "differential"), dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), weighted=FALSE) {
		type <- match.arg(type)
		dose <- match.arg(dose)
		dose.units <- match.arg(dose.units)
		volume <- match.arg(volume)
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return(x)
		}
		x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
		structure.name <- paste("mean(", paste(names(x), collapse=", ", sep=""), ")", sep="")
		structure.volumes <- as.numeric(lapply(x, slot, "structure.volume"))
		structure.means <- as.numeric(lapply(x, slot, "dose.mean"))
		if (weighted) {
			dose.mean <- sum(structure.means * structure.volumes, na.rm=TRUE)/sum(structure.volumes, na.rm=TRUE)
		}
		else {
			dose.mean <- mean(structure.means, na.rm=TRUE)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		dose.rx <- max(100 * unlist(lapply(x,slot,"dose.rx")) / unlist(lapply(x,slot,"rx.isodose")), na.rm=TRUE)
		rx.isodose <- 100
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		if (weighted) {
			volumes.new <- apply(t(t(volume.matrix)*structure.volumes), 1, sum, na.rm=TRUE)/sum(structure.volumes)
		}
		else {
			volumes.new <- apply(volume.matrix, 1, mean, na.rm=TRUE)
		}
		new("DVH", type=type, dose.type=dose, volume.type=volume, structure.name=structure.name, structure.volume=mean(structure.volumes), dose.min=dose.min, dose.rx=dose.rx, rx.isodose=rx.isodose, dose.max=dose.max, dose.mean=dose.mean, doses=doses.new, dose.units=dose.units, volumes=volumes.new)
	}
)


setMethod("mean", "DVH",
	function (x, na.rm=TRUE) {
		if (is.na(x@dose.mean)) {
			if (x@type == "cumulative") {
				x <- convert.DVH(x, type="differential", dose=x@dose.type, volume="absolute", dose.units=x@dose.units)
			}
			return(sum(x@doses*x@volumes, na.rm=na.rm)/sum(x@volumes, na.rm=na.rm))
		}
		else {
			return(x@dose.mean)
		}
	}
)


setGeneric("median",
	median
)

setMethod("median", "DVH.list",
	function (x, na.rm=TRUE) {
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return(x)
		}
		identical.format <- TRUE
		type.format <- unique(unlist(lapply(x, slot, "type")))
		if (length(type.format) > 1) {
			warning("DVH 'type' parameter differs among DVH data (will reinterpet all data as cumulative DVH)")
			type.format <- "cumulative"
			identical.format <- FALSE
		}
		dose.format <- unique(unlist(lapply(x, slot, "dose.type")))
		if (length(dose.format) > 1) {
			warning("'dose.type' parameter differs among DVH data (will reinterpret all data relative dose)")
			dose.format <- "relative"
			identical.format <- FALSE
		}
		volume.format <- unique(unlist(lapply(x, slot, "volume.type")))
		if (length(volume.format) > 1) {
			warning("'volume.type' parameter differs among DVH data (will reinterpret all data as relative volume)")
			volume.format <- "relative"
			identical.format <- FALSE
		}
		units.format <- unique(unlist(lapply(x, slot, "dose.units")))
		if (length(units.format) > 1) {
			warning("'dose.units' parameter differs among DVH data (will reinterpret all doses as cGy)")
			units.format <- "cGy"
			identical.format <- FALSE
		}
		if (!identical.format) {
			x <- convert.DVH(x, type=type.format, dose=dose.format, volume=volume.format, dose.units=units.format)
		}
		structure.name <- paste("median(", paste(names(x), collapse=", ", sep=""), ")", sep="")
		structure.volume <- median(as.numeric(lapply(x, slot, "structure.volume")), na.rm=na.rm)
		structure.means <- as.numeric(lapply(x, slot, "dose.mean"))
		dose.mean <- median(structure.means, na.rm=na.rm)
		dose.units <- slot(x[[1]], "dose.units")
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		dose.rx <- max(100 * unlist(lapply(x,slot,"dose.rx")) / unlist(lapply(x,slot,"rx.isodose")), na.rm=TRUE)
		rx.isodose <- 100
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		volumes.new <- apply(volume.matrix, 1, median, na.rm=na.rm)
		new("DVH", type=x[[1]]$type, dose.type=x[[1]]$dose.type, volume.type=x[[1]]$volume.type, structure.name=structure.name, structure.volume=structure.volume, dose.rx=dose.rx, rx.isodose=rx.isodose, dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, doses=doses.new, dose.units=dose.units, volumes=volumes.new)
	}
)


setGeneric("mad", 
	mad
)

setMethod("mad", "DVH.list",
	function (x) {
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return()
		}
		identical.format <- TRUE
		type.format <- unique(unlist(lapply(x, slot, "type")))
		if (length(type.format) > 1) {
			warning("DVH 'type' parameter differs among DVH data (will reinterpet all data as cumulative DVH)")
			type.format <- "cumulative"
			identical.format <- FALSE
		}
		dose.format <- unique(unlist(lapply(x, slot, "dose.type")))
		if (length(dose.format) > 1) {
			warning("'dose.type' parameter differs among DVH data (will reinterpret all data relative dose)")
			dose.format <- "relative"
			identical.format <- FALSE
		}
		volume.format <- unique(unlist(lapply(x, slot, "volume.type")))
		if (length(volume.format) > 1) {
			warning("'volume.type' parameter differs among DVH data (will reinterpret all data as relative volume)")
			volume.format <- "relative"
			identical.format <- FALSE
		}
		units.format <- unique(unlist(lapply(x, slot, "dose.units")))
		if (length(units.format) > 1) {
			warning("'dose.units' parameter differs among DVH data (will reinterpret all doses as cGy)")
			units.format <- "cGy"
			identical.format <- FALSE
		}
		if (!identical.format) {
			x <- convert.DVH(x, type=type.format, dose=dose.format, volume=volume.format, dose.units=units.format)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, mad=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, mad=apply(volume.matrix, 1, mad)))
	}
)


setGeneric("quantile", 
	quantile
)

setMethod("quantile", "DVH.list",
	function (x, type=7, ...) {
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return()
		}
		identical.format <- TRUE
		type.format <- unique(unlist(lapply(x, slot, "type")))
		if (length(type.format) > 1) {
			warning("DVH 'type' parameter differs among DVH data (will reinterpet all data as cumulative DVH)")
			type.format <- "cumulative"
			identical.format <- FALSE
		}
		dose.format <- unique(unlist(lapply(x, slot, "dose.type")))
		if (length(dose.format) > 1) {
			warning("'dose.type' parameter differs among DVH data (will reinterpret all data relative dose)")
			dose.format <- "relative"
			identical.format <- FALSE
		}
		volume.format <- unique(unlist(lapply(x, slot, "volume.type")))
		if (length(volume.format) > 1) {
			warning("'volume.type' parameter differs among DVH data (will reinterpret all data as relative volume)")
			volume.format <- "relative"
			identical.format <- FALSE
		}
		units.format <- unique(unlist(lapply(x, slot, "dose.units")))
		if (length(units.format) > 1) {
			warning("'dose.units' parameter differs among DVH data (will reinterpret all doses as cGy)")
			units.format <- "cGy"
			identical.format <- FALSE
		}
		if (!identical.format) {
			x <- convert.DVH(x, type=type.format, dose=dose.format, volume=volume.format, dose.units=units.format)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, quantiles=apply(volume.matrix, 1, quantile, type=type, ...)))
	}
)

setGeneric("var", 
	var
)

setMethod("var", "DVH.list",
	function (x, na.rm=TRUE) {
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return()
		}
		identical.format <- TRUE
		type.format <- unique(unlist(lapply(x, slot, "type")))
		if (length(type.format) > 1) {
			warning("DVH 'type' parameter differs among DVH data (will reinterpet all data as cumulative DVH)")
			type.format <- "cumulative"
			identical.format <- FALSE
		}
		dose.format <- unique(unlist(lapply(x, slot, "dose.type")))
		if (length(dose.format) > 1) {
			warning("'dose.type' parameter differs among DVH data (will reinterpret all data relative dose)")
			dose.format <- "relative"
			identical.format <- FALSE
		}
		volume.format <- unique(unlist(lapply(x, slot, "volume.type")))
		if (length(volume.format) > 1) {
			warning("'volume.type' parameter differs among DVH data (will reinterpret all data as relative volume)")
			volume.format <- "relative"
			identical.format <- FALSE
		}
		units.format <- unique(unlist(lapply(x, slot, "dose.units")))
		if (length(units.format) > 1) {
			warning("'dose.units' parameter differs among DVH data (will reinterpret all doses as cGy)")
			units.format <- "cGy"
			identical.format <- FALSE
		}
		if (!identical.format) {
			x <- convert.DVH(x, type=type.format, dose=dose.format, volume=volume.format, dose.units=units.format)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, var=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- pmax(0, approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y, na.rm=TRUE)
		}
		return(list(dose=doses.new, var=apply(volume.matrix, 1, var, na.rm=na.rm)))
	}
)

setGeneric("sd", 
	sd
)

setMethod("sd", "DVH.list",
	function (x, na.rm=TRUE) {
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N < 1) {
			return()
		}
		identical.format <- TRUE
		type.format <- unique(unlist(lapply(x, slot, "type")))
		if (length(type.format) > 1) {
			warning("DVH 'type' parameter differs among DVH data (will reinterpet all data as cumulative DVH)")
			type.format <- "cumulative"
			identical.format <- FALSE
		}
		dose.format <- unique(unlist(lapply(x, slot, "dose.type")))
		if (length(dose.format) > 1) {
			warning("'dose.type' parameter differs among DVH data (will reinterpret all data relative dose)")
			dose.format <- "relative"
			identical.format <- FALSE
		}
		volume.format <- unique(unlist(lapply(x, slot, "volume.type")))
		if (length(volume.format) > 1) {
			warning("'volume.type' parameter differs among DVH data (will reinterpret all data as relative volume)")
			volume.format <- "relative"
			identical.format <- FALSE
		}
		units.format <- unique(unlist(lapply(x, slot, "dose.units")))
		if (length(units.format) > 1) {
			warning("'dose.units' parameter differs among DVH data (will reinterpret all doses as cGy)")
			units.format <- "cGy"
			identical.format <- FALSE
		}
		if (!identical.format) {
			x <- convert.DVH(x, type=type.format, dose=dose.format, volume=volume.format, dose.units=units.format)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, sd=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, sd=apply(volume.matrix, 1, sd, na.rm=na.rm)))
	}
)

setMethod("max", "DVH.list",
	function(x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(max(as.numeric(lapply(x, max, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("max", "DVH",
	function (x, na.rm=TRUE) {
		if (is.na(x@dose.max)) {
			return(max(x@doses, na.rm=na.rm))
		}
		else {
			return(x@dose.max)
		}
	}
)

setMethod("min", "DVH.list",
	function(x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(min(as.numeric(lapply(x, min, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("min", "DVH",
	function (x, na.rm=TRUE) {
		if (is.na(x@dose.min)) {
			return(min(x@doses, na.rm=na.rm))
		}
		else {
			return(x@dose.min)
		}
	}
)


setMethod("range", "DVH.list",
	function (x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(range(unlist(lapply(x, range, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("range", "DVH",
	function (x, ..., na.rm=TRUE) {
		if (is.na(x@dose.min) | is.na(x@dose.max)) {
			return(range(x@doses, na.rm=na.rm))
		}
		else {
			return(c(x@dose.min, x@dose.max))
		}
	}
)


setMethod("sum", "DVH.list",
	function (x, ..., type=c("cumulative", "differential"), dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), na.rm = TRUE) {
		x <- c(x, ...)
		x <- x[!unlist(lapply(x, is.empty))]
		N <- length(x)
		if (N <= 1) {
			return(x)
		}
		type <- match.arg(type)
		dose <- match.arg(dose)
		dose.units <- match.arg(dose.units)
		volume <- match.arg(volume)
		x <- new("DVH.list", lapply(x, convert.DVH, type="differential", dose="absolute", dose.units=dose.units, volume="absolute"))
		structure.name <- paste("sum(", paste(names(x), collapse=", ", sep=""), ")", sep="")
		structure.volumes <- as.numeric(lapply(x, slot, "structure.volume"))
		volume.sum <- sum(structure.volumes, na.rm=TRUE)
		dose.max <- max(x, na.rm=na.rm)
		dose.min <- min(x, na.rm=na.rm)
		dose.mean <- as.numeric(lapply(x, slot, "dose.mean"))
		dose.mean <- sum(dose.mean * structure.volumes, na.rm=TRUE)/volume.sum
		dose.rx <- max(100 * unlist(lapply(x,slot,"dose.rx")) / unlist(lapply(x,slot,"rx.isodose")), na.rm=TRUE)
		rx.isodose <- 100
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		doses.new <- unique(unlist(lapply(x, slot, "doses")))
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			if (class(x[[i]]) == "zDVH") {
				x[[i]] <- as(x[[i]], "DVH")	
			}
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		volumes.new <- apply(volume.matrix, 1, sum, na.rm=na.rm)
		x <- new("DVH", type="differential", dose.type="absolute", volume.type="absolute", structure.name=structure.name, structure.volume=volume.sum, dose.min=dose.min, dose.rx=dose.rx, rx.isodose=rx.isodose, dose.max=dose.max, dose.mean=dose.mean, doses=doses.new, dose.units=dose.units, volumes=volumes.new)
		return(convert.DVH(x, type=type, dose=dose, dose.units=dose.units, volume=volume))
	}
)