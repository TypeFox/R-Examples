### Written by Romain Francois <francoisromain@free.fr>
### TODO: allow also to indicate imported packages, may be?
pkgLoaded <- function ()
{
	s <- grep("^package:", search(), value = TRUE)
	return(sub("^package:", "", s))
}

### Written by Romain Francois <francoisromain@free.fr>
pkgInstalled <- function (pattern = NULL, ...)
{
	## Warning: installed.packages() can be very slow!
	ip <- installed.packages(fields = "Title", ...)
	if (!is.null(pattern)) {
		keep <- suppressWarnings(union( 
			grep(pattern, ip [, "Package"], ignore.case = TRUE), 
			grep(pattern, ip [, "Title"], ignore.case = TRUE)))
		ip <- ip[keep, , drop = FALSE]
	}
	lp <- pkgLoaded() 
	def <- c(getOption("defaultPackages"), "base")
	ip <- cbind(ip, 
		"Loaded" = ip[, "Package"] %in% lp, 
		"Default" = ip[, "Package"] %in% def)
	return(ip) 
}

### Written by Romain Francois <francoisromain@free.fr>
pkgDesc <- function (pkg, lib.loc = NULL, fields = NULL, encoding = "")
{
    retval <- list()
    if (!is.null(fields)) {
        fields <- as.character(fields)
        retval[fields] <- NA
    }
    pkgpath <- ""
    if (is.null(lib.loc)) {
        if (pkg == "base") 
            pkgpath <- file.path(.Library, "base")
        else if ((envname <- paste("package:", pkg, sep = "")) %in% search()) {
            pkgpath <- attr(as.environment(envname), "path")
            if (is.null(pkgpath)) pkgpath <- ""
        }
    }
    if (pkgpath == "") {
        libs <- if (is.null(lib.loc)) .libPaths() else lib.loc
        for (lib in libs) {
			if (file.access(file.path(lib, pkg), 5) == 0) {
				pkgpath <- file.path(lib, pkg)
				break
			}
		}
    }
    if (pkgpath == "") {
        pkgpath <- system.file(package = pkg, lib.loc = lib.loc)
        if (pkgpath == "") {
            warning(gettextf("no package '%s' was found", pkg), domain = NA)
            return(NA)
        }
    }
    file <- file.path(pkgpath, "DESCRIPTION")
    return(readLines(file))
}

### Written by Romain Francois <francoisromain@free.fr>
pkgWebDesc <- function (pkg, repos = getOption("repos"))
{
	temp <- tempfile()
	on.exit(unlink(temp))
	txt <- suppressWarnings(try({
		download.file(sprintf("%s/web/packages/%s/DESCRIPTION", repos, pkg),
			destfile = temp, quiet = TRUE)
		readLines(temp)
	}, silent = TRUE))
	if (inherits(txt, "try-error")) txt <- ""
	return(txt)
}
