make_EPSG <- function(file) {
        if (missing(file)) {
            tf <- tempfile()
            n <- .Call("PROJcopyEPSG", tf, PACKAGE="rgdal")
            if (n > 0) file <- tf
            else stop("Error opening epsg file")
       } else {
# bug fix 20150914 Alain Benard
           if (file.exists(file)) {
               n <- length(readLines(file))
           } else {
               stop("Error opening", file)
           }
       }
#	epsg <- file(system.file("proj/epsg", package="rgdal"), open="r")
# report Michael Sumner 071017
	epsg <- file(file, open="r")
	note <- character(n)
	prj4 <- character(n)
	code <- integer(n)
	for (i in seq(along=code)) {
		res <- strsplit(paste(readLines(epsg, n=2), collapse=""), 
			"[<>]")
		note[i] <- res[[1]][1]
		prj4[i] <- res[[1]][3]
		code[i] <- as.integer(res[[1]][2])
		if (is.na(code[i])) res <- readLines(epsg, n=1)
	}
	close(epsg)
        prj4 <- scan(con <- textConnection(prj4), what=character(0), sep="\n",
                 strip.white=TRUE, quiet=TRUE)
        close(con)
#	prj4 <- scan(textConnection(prj4), what=character(0), sep="\n", 
#		strip.white=TRUE, quiet=TRUE)
	EPSG <- data.frame(code=I(code), note=I(note), prj4=I(prj4))
	EPSG <- EPSG[-which(is.na(EPSG$note)),]
	EPSG
}

