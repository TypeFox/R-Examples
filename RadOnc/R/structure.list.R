
setMethod("as.list", "structure.list",
	function(x, ...) {
		return(attr(x,"structures"))
	}
)

setAs("list", "structure.list", 
	function(from) {
		struct.list.combined <- new("structure.list")
		lapply(from, function (struct.list) {
			struct.list.combined <<- c(struct.list.combined, struct.list)
		})
		return(struct.list.combined)
	}
)

setMethod("lapply", "structure.list",
	function (X, FUN, ...) {
    	X <- as.list(X)
    	.Internal(lapply(X, FUN))
	}
)


setMethod("length", "structure.list",
	function (x) {
		return(length(attr(x,"structures")))
	}
)


setMethod("[", "structure.list",
	function (x, i, ...) {
		if (missing(i) || (length(i) < 1) || is.na(i)) {
			return(new("structure.list"))
		}
		if (all(is.logical(i))) {
			x <- attr(x,"structures")
			return(new("structure.list", x[i]))
		}
		if (suppressWarnings(all(!is.na(as.numeric(i))))) {
			x <- attr(x,"structures")
			return(new("structure.list", x[as.numeric(i)]))
		}
		if (length(i) == 1) {
			names.x <- names(x)
			x <- attr(x,"structures")
			if (grepl("(\\*|\\^|\\$|\\?|\\+|[[]|[{]|\\|)", i)) {
				return(new("structure.list", x[grep(i, names.x)]))
			}
			else if (is.character(i)) {
				return(new("structure.list", x[which(names.x == i)]))
			}
			else if (is.logical(i)) {
				return(new("structure.list", x[i]))
			}			
			else if (suppressWarnings(!is.na(as.numeric(i)))) {
				return(new("structure.list", x[as.numeric(i)]))
			}
			else {
				return(new("structure.list", x[i]))
			}			
		}
		return(c(x[i[1]], x[i[2:length(i)]]))
	}
)


setMethod("$", "structure.list",
	function (x, name) {
		name <- unlist(strsplit(name, ","))
		return(lapply(x, function (struct) { struct[name] }))		
	}
)


setMethod("[[", "structure.list",
	function (x, i, exact=TRUE) {
		x <- attr(x,"structures")
		return(x[[i]])
	}
)


setMethod("[[<-", "structure.list",
	function (x, i, value) {
		x <- attr(x,"structures")
		if (class(value) == "structure3D") {
			x[[i]] <- value
		}
		else {
			stop("'value' must be an object of class 'structure3D'")
		}
		return(new("structure.list", x))
	}
)

setMethod("c", "structure.list",
	function (x, ..., recursive = FALSE) {
		return(new("structure.list", c(as.list(x), as.list(c(... , recursive=FALSE)))))
	}
)


setMethod("rev", "structure.list",
	function (x) {
		if (length(x) <= 1) {
			return(x)
		}
		else {
			return(x[length(x):1])
		}
	}
)


setMethod("print", "structure.list",
	function (x, ...) {
		print(paste("List containing ", length(x), " structure3D objects (", paste(names(x), collapse=", ", sep=""), ")", sep=""))
	}
)


setMethod("show", "structure.list",
	function (object) {
		print(object)
	}
)


setMethod("names", "structure.list",
	function (x) {
		return(as.character(unlist(lapply(x, names))))
	}
)


setMethod("names<-", "structure.list",
 	function (x, value) {
		if (length(x) != length(value)) {
			stop(paste("'names' attribute [", length(value), "] must be the same length as the structure3D list [", length(x), "]", sep=""))
		}
		struct.list <- new("structure.list", mapply(function(struct, name) {
				struct$name <- name
				return(struct)
			},
			x, value
		))
		names(attr(struct.list,"structures")) <- value
		return(struct.list)
  	}
)

setMethod("range", "structure.list",
	function (x, ..., na.rm=TRUE) {
		if (length(x) < 1) {
			warning("Cannot calculate range of an empty structure list")
			return(matrix(NA, nrow=2, ncol=3, dimnames=list(c("min", "max"), c("x", "y", "z"))))
		}
		ranges <- lapply(x, range)
		range <- matrix(rep(c(Inf, -Inf), 3), nrow=2, ncol=3, dimnames=list(c("min", "max"), c("x", "y", "z")))
		for (i in 1:length(ranges)) {
			range[1, ] <- pmin(range[1, ], ranges[[i]][1, ], na.rm=na.rm)
			range[2, ] <- pmax(range[2, ], ranges[[i]][2, ], na.rm=na.rm)
		}
		return(range)
	}
)


setMethod("plot", c("structure.list", "missing"),
	function(x, y, col="gray", alpha=1, ...) {
		open3d()
		if (length(x) != length(col)) {
			col <- rep(col[1], length(x))
		}
		if (length(x) != length(alpha)) {
			alpha <- rep(alpha[1], length(x))
		}
		value <- mapply(
			function(x, col, alpha) {
				if (dim(x$triangles)[2] >= 1) {
					triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
				}
				else {
					points3d(x$vertices, col=col, alpha=alpha)
				}
			},
			x,
			col,
			alpha
		)
		
	}
)

setMethod("plot", c("structure.list", "ANY"),
	function(x, y, col="gray", alpha=1, ...) {
		open3d()
		if (length(x) != length(col)) {
			col <- rep(col[1], length(x))
		}
		if (length(x) != length(alpha)) {
			alpha <- rep(alpha[1], length(x))
		}
		value <- mapply(
			function(x, col, alpha) {
				if (dim(x$triangles)[2] >= 1) {
					triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
				}
				else {
					points3d(x$vertices, col=col, alpha=alpha)
				}
			},
			x,
			col,
			alpha
		)
		
	}
)
