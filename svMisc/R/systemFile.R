systemFile <- function (..., exec = FALSE, package = NULL, lib.loc = NULL)
{
	## First look if exec is TRUE
	if (isTRUE(exec)) {
		res <- Sys.which(as.character(unlist(list(...))))
		if (length(res) == 1) res <- as.character(res)
	} else if (!is.null(package)) {
		## A file in a package
		res <- system.file(..., package = package, lib.loc = lib.loc)
		## Check that this is a directory
		if (!file_test("-f", res)) res <- ""
	} else {
		## Look if this file exists and is a file
		file <- as.character(unlist(list(...)))
		file <- file.path(file)
		if (file_test("-f", file)) res <- normalizePath(file) else res <- ""
	}
	res
}

systemDir <- function (..., exec = FALSE, package = NULL, lib.loc = NULL)
{
	## First look if exec is TRUE
	if (isTRUE(exec)) {
		files <- Sys.which(as.character(unlist(list(...))))
		## Note: Sys.which() does not always return "" for items not found!
		res <- dirname(files)
		res[res == "."] <- ""
		if (length(res) > 1) names(res) <- names(files)
	} else if (!is.null(package)) {
		## A directory in a package
		res <- system.file(..., package = package, lib.loc = lib.loc)
		## Check that this is a directory
		if (!file_test("-d", res)) res <- ""
	} else {
		## A predefined directory
		which <- as.character(unlist(list(...)))
		
		## This is a specific directory
		getDir <- function (which = c("temp", "sysTemp", "user", "home", "bin",
			"doc", "etc", "share")) {
			which = match.arg(which)
			res <- switch(which,
				"temp" = tempdir(),
				"sysTemp" = if (!isWin() && file_test("-d", "/tmp")) "/tmp" else
					dirname(tempdir()),
				"user" = tools::file_path_as_absolute("~"),
				"home" = R.home("home"),
				"bin" = R.home("bin"),
				"doc" = R.home("doc"),
				"etc" = R.home("etc"),
				"share" = R.home("share"))
			return(res)
		}
		if (is.null(which) || length(which) == 0) return(character(0)) else {
			res <- character(length(which))
			if (length(which) > 1) names(res) <- which
			for (i in seq_along(which))
				res[i] <- getDir(which[i])
		}
	}
	res
}
