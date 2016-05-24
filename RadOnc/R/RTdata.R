setMethod("print", "RTdata",
	function (x, ...) {
		contents <- c()
		if (any(dim(slot(x, "CT")) > 0)) {
			contents <- c(contents, paste("CT image (", paste(dim(slot(x, "CT")), collapse="x", sep=""), ")", sep=""))
		}
		if (any(dim(slot(x, "dose")) > 0)) {
			contents <- c(contents, paste("dose grid (", paste(dim(slot(x, "dose")), collapse="x", sep=""), ")", sep=""))
		}
		if (length(slot(x, "structures")) > 0) {
			contents <- c(contents, paste(length(slot(x, "structures")), "structure(s)"))
		}
		if (length(contents) == 0) {
			print("Empty RT dataset")
		}
		else {
			print(paste("RT data '", slot(x, "name"), "' containing ", paste(contents, collapse=", ", sep=""), sep=""))
		}
	}
)


setMethod("show", "RTdata",
	function (object) {
		print(object)
	}
)


setMethod("$", "RTdata",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'RTdata'")
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("names", "RTdata",
	function (x) {
		return(x$name)
	}
)


setMethod("names<-", "RTdata",
 	function (x, value) {
 		x$name <- value
 		return(x)
 	}
)


setMethod("$<-", "RTdata",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'RTdata'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)
