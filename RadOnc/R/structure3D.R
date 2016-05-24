setMethod("$", "structure3D",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'structure3D'")
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("names", "structure3D",
	function (x) {
		return(x$name)
	}
)


setMethod("names<-", "structure3D",
 	function (x, value) {
 		x$name <- value
 		return(x)
 	}
)


setMethod("$<-", "structure3D",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'structure3D'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)

setMethod("c", "structure3D",
	function (x, ..., recursive = FALSE) {
		return(c(as(x, "structure.list"), ..., recursive=FALSE))
	}
)


setMethod("range", "structure3D",
	function (x, ..., na.rm=TRUE) {
		if (dim(x$vertices)[1] >= 1) {
			range <- apply(x$vertices, 2, range, na.rm=na.rm)
		}
		else {
			range <- matrix(NA, nrow=2, ncol=3)
		}
		dimnames(range) <- list(c("min", "max"), c("x", "y", "z"))
		return(range)
	}
)


setMethod("plot", c("structure3D", "missing"),
	function(x, col="gray", alpha=1, ...) {
		open3d()
		if (dim(x$triangles)[2] >= 1) {
			triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
		}
		else {
			points3d(x$vertices, col=col, alpha=alpha)
		}
	}
)

setMethod("plot", c("structure3D", "ANY"),
	function(x, y, col="gray", alpha=1, ...) {
		open3d()
		if (dim(x$triangles)[2] >= 1) {
			triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
		}
		else {
			points3d(x$vertices, col=col, alpha=alpha)
		}
	}
)

setAs("structure3D", "structure.list", 
	function(from) {
		return(new("structure.list", structures=from))
	}
)

setMethod("print", "structure3D",
	function (x, ...) {
		print(paste("Structure (", names(x), ") defined by ", dim(x$vertices)[1], " points in ", length(x$closed.polys), " axial slices", sep=""))
	}
)


setMethod("show", "structure3D",
	function (object) {
		print(object)
	}
)

setMethod("dim", "structure3D",
	function (x) {
		return(c(dim(attr(x,"vertices"))[1], length(attr(x,"closed.polys"))))
	}
)
