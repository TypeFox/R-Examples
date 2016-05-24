
setMethod("$", "zDVH",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)

setMethod("names", "zDVH",
	function (x) {
		return(x$structure.name)
	}
)

setMethod("names<-", "zDVH",
 	function (x, value) {
 		x$structure.name <- value
 		return(x)
 	}
)

setMethod("$<-", "zDVH",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'zDVH'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)

setMethod("[", "zDVH",
	function (x, i, ...) {
		if (!validObject(x)) {
			stop("not a valid object of 'zDVH'")
		}
		x <- as(x, "DVH")
		return(x[i])
	}
)


setMethod("c", "zDVH",
	function (x, ..., recursive = FALSE) {
		return(c(as(x, "DVH.list"), ..., recursive=recursive))
	}
)

setMethod("sum", "DVH",
	function (x, ..., na.rm = TRUE) {
		return(sum(as(x, "DVH.list"), ..., na.rm = na.rm))
	}
)

setMethod("print", "zDVH",
	function (x, ...) {
		if (x@dose.type == "relative") {
			dose.type <- "%"
			dose.min <- x@dose.min * x@rx.isodose / x@dose.rx
			dose.max <- x@dose.max * x@rx.isodose / x@dose.rx
		}
		else {
			dose.type <- x@dose.units
			dose.min <- x@dose.min
			dose.max <- x@dose.max
		}
		print(paste("Structure: ", x@structure.name, " (", sprintf("%.*f", 1, x@structure.volume), "cc), Dose: ", sprintf("%.*f", 2, dose.min), "-", sprintf("%.*f", 2, dose.max), dose.type, " (", x@dose.rx, x@dose.units, " prescribed", if (x@rx.isodose != 100) {paste(" to ", x@rx.isodose, "% isodose line", sep="")}, "), DVH: ", x@type, ", Volume: ", x@volume.type, ", Axial segments: ", dim(x@volumes)[2], sep=""))
	}
)

setMethod("show", "zDVH",
	function (object) {
		print(object)
	}
)
