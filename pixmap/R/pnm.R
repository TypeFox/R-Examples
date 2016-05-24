read.pnm <- function(file, ...)
{
    fsz <- file.info(file)$size
    con <- file(file, open="rb")

    pnmhead <- read.pnmhead(con)
    retval <- read.pnmdata(con, pnmhead, ...)

    if (fsz != seek(con))
	warning("Possible reading error: file size ", fsz,
		" bytes, but ", seek(con), " bytes read")
    close(con)
    retval
}

read.pnmdata <- function(con, pnmhead, ...)
{
    ds <- pnmhead$datastart
    seek(con, ds)
    type <- pnmhead$type
    nl <- ifelse(type == "ppm", 3, 1)
    nc <- pnmhead$nc
    nr <- pnmhead$nr
    ncells <- nl*nc*nr
    if (pnmhead$ascii) {
	xx <- scan(con, integer(0), n=ncells)
    }
    else {
	if (type == "pbm") { ## black & white, i.e. pixel = bit
	    BytesPerRow <- ceiling(nc/8)
	    bxx <- readBin(con, "integer",
			   n=nr*BytesPerRow, size=1, signed=FALSE)

	    as.integer.bytes <- function (x) {
		## unpacks bytes in 0:255 into {0,1} integers
		n <- length(x <- as.integer(x))
		if (any(x < 0) || any(x > 255))
		    stop("Not an unsigned byte (value outside 0:255)")
		ans <- matrix(integer(8 * n), 8, n)
		two <- as.integer(2)
		for (i in 8:1) {
		    ans[i,] <- x %% two
		    x <- x %/% two
		}
		ans
	    }

	    xx <- as.integer.bytes(bxx)
	    ncb <- BytesPerRow*8
	    xx <- 1 - array(xx, c(nl, ncb, nr))[,1:nc,]
	}
	else {
	    xx <- readBin(con, "integer",
			  n=ncells, size=1, signed=FALSE)
	}
    }

    res <- array(xx, dim = c(nl, nc, nr)) / pnmhead$maxval

    if(nl==1) { ## non-RGB:
	##FIXME(MM): use  "indexed" for B&W
	z = pixmapGrey(t(res[1,,]), ...)
    }
    else{
        z = pixmapRGB(0, ncol=dim(res)[2], nrow=dim(res)[3], ...)
        z@red = t(res[1,,])
        z@green = t(res[2,,])
        z@blue = t(res[3,,])
    }
    z
}



write.pnm <- function(object, file=NULL, forceplain=FALSE,
                      type=NULL, maxval=255)
{
    if(!is(object, "pixmap"))
        stop("Can only write pixmap objects")

    if(is.null(type))
        type <- if(is(object, "pixmapGrey")) "pgm" else "ppm"
    else
        type <- match.arg(type, c("pbm", "pgm", "ppm"))
    do <- object@size

    switch(type,
           "pbm" = {
               object <- as(object, "pixmapGrey")
               object <- t(object@grey < 0.5)
               storage.mode(object) <- "integer"
               code <- 4
               forceplain <- TRUE
           },
           "pgm" = {
               object <- as(object, "pixmapGrey")
               object <- t(round(object@grey*maxval, 0))
               storage.mode(object) <- "integer"
               code <- 5
           },
           "ppm" = {
               object <- as(object, "pixmapRGB")
               object1 <- array(0, dim=c(3, do[2], do[1]))
               object1[1,,] <- t(object@red)
               object1[2,,] <- t(object@green)
               object1[3,,] <- t(object@blue)
               object <- object1
               object <- round(object*maxval, 0)
               storage.mode(object) <- "integer"
               code <- 6
           })

    if (is.null(file)) file <- paste("Rimage.", type, sep="")
    comment <- "# R write.pnm output"
    if(forceplain) {
        con <- file(file, open="w")
        code <- code - 3
        cat("P", code, "\n", file=con, sep="")
        cat(comment, "\n", file=con, sep="")
        cat(do[2], " ", do[1], "\n", file=con, sep="")
        if (type != "pbm")
            cat(maxval, "\n", file=con, sep="")

        write(object, ncolumns=3, file=con)
    }
    else {
        con <- file(file, open="wb")
        writeChar(paste("P", code, "\n", sep=""), con=con, eos=NULL)
        writeChar(paste(comment, "\n", sep=""), con=con, eos=NULL)
        writeChar(paste(do[2], " ", do[1], "\n", sep=""),
                  con=con, eos=NULL)
        if (type != "pbm") writeChar(paste(maxval, "\n", sep=""),
              con=con, eos=NULL)
        writeBin(as.integer(as.vector(object)), con, size=1)
    }
    close(con)
}

