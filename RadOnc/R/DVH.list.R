
setMethod("as.list", "DVH.list",
	function(x, ...) {
		return(attr(x,"structures"))
	}
)

setAs("DVH", "DVH.list", 
	function(from) {
		return(new("DVH.list", structures=from))
	}
)

setAs("zDVH", "DVH.list", 
	function(from) {
		return(new("DVH.list", structures=from))
	}
)

setAs("structure.list", "DVH.list", 
	function(from) {
		return(new("DVH.list", lapply(from, function(struct) {return(struct$DVH)})))
	}
)

setAs("list", "DVH.list", 
	function(from) {
		DVH.list.combined <- new("DVH.list")
		lapply(from, function (DVH.list) {
			DVH.list.combined <<- c(DVH.list.combined, DVH.list)
		})
		return(DVH.list.combined)
	}
)

setMethod("lapply", "DVH.list",
	function (X, FUN, ...) {
    	X <- as.list(X)
    	.Internal(lapply(X, FUN))
	}
)


setMethod("length", "DVH.list",
	function (x) {
		return(length(attr(x,"structures")))
	}
)


setMethod("[", "DVH.list",
	function (x, i, ...) {
		if (missing(i) || (length(i) < 1) || is.na(i)) {
			return(new("DVH.list"))
		}
		if (all(is.logical(i))) {
			x <- attr(x,"structures")
			return(new("DVH.list", x[i]))
		}
		if (suppressWarnings(all(!is.na(as.numeric(i))))) {
			x <- attr(x,"structures")
			return(new("DVH.list", x[as.numeric(i)]))
		}
		if (length(i) == 1) {
			x <- attr(x,"structures")
			if (grepl("(\\*|\\^|\\$|\\?|\\+|[[]|[{]|\\|)", i)) {
				return(new("DVH.list", x[grep(i, unlist(lapply(x, names)))]))
			}
			else if (is.character(i)) {
				return(new("DVH.list", x[which(unlist(lapply(x, names)) == i)]))					
			}
			else if (is.logical(i)) {
				return(new("DVH.list", x[i]))
			}
			else if (suppressWarnings(!is.na(as.numeric(i)))) {
				return(new("DVH.list", x[i]))
			}
			else {
				return(new("DVH.list", x[i]))
			}			
		}
		return(c(x[i[1]], x[i[2:length(i)]]))
	}
)

setMethod("$", "DVH.list",
	function (x, name) {
		name <- unlist(strsplit(name, ","))
		return(lapply(x, function (DVH) { DVH[name] }))		
	}
)


setMethod("[[", "DVH.list",
	function (x, i, exact=TRUE) {
		x <- attr(x,"structures")
		return(x[[i]])
	}
)

setMethod("[[<-", "DVH.list",
	function (x, i, value) {
		x <- attr(x,"structures")
		if (class(value) %in% c("DVH", "zDVH")) {
			x[[i]] <- value
		}
		else {
			stop("'value' must be an object of class 'DVH' or 'zDVH'")
		}
		return(new("DVH.list", x))
	}
)

setMethod("c", "DVH.list",
	function (x, ..., recursive = FALSE) {
		return(new("DVH.list", c(as.list(x), as.list(c(... , recursive=FALSE)))))
	}
)

setMethod("rev", "DVH.list",
	function (x) {
		if (length(x) <= 1) {
			return(x)
		}
		else {
			return(x[length(x):1])
		}
	}
)

setMethod("print", "DVH.list",
	function (x, ...) {
		print(paste("List containing ", length(x), " DVH objects (", paste(names(x), collapse=", ", sep=""), ")", sep=""))
	}
)

setMethod("show", "DVH.list",
	function (object) {
		print(object)
	}
)

setMethod("names", "DVH.list",
	function (x) {
		return(as.character(unlist(lapply(x, names))))
	}
)

setMethod("names<-", "DVH.list",
 	function (x, value) {
		if (length(x) != length(value)) {
			stop(paste("'names' attribute [", length(value), "] must be the same length as the DVH list [", length(x), "]", sep=""))
		}
		DVHlist <- new("DVH.list", mapply(function(DVH, name) {
				DVH$structure.name <- name
				return(DVH)
			},
			x, value
		))
		names(attr(DVHlist,"structures")) <- value		
		return(DVHlist)
  	}
)
