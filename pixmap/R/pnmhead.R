read.pnmhead <- function(con) {

    seek(con, 0)

    pm.getc <- function(con) {
	ch <- readChar(con, nchars=1)
	if (ch == "#") {
	    ch <- readChar(con, nchars=1)
	    while (ch != '\n' && ch != '\r') {
		ch <- readChar(con, nchars=1)
	    }
	}
	ch	
    }

    pm.getuint <- function(con) {
	ch <- pm.getc(con)
	while (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r') {
	    ch <- pm.getc(con)
	}
	if (ch < '0' || ch > '9')
	    stop("junk in file where an unsigned integer should be")
	i <- 0
	while (ch >= '0' && ch <= '9') {
	    digitVal <- as.integer(ch)
	    i <- i * 10 + digitVal
	    ch <- pm.getc(con)
	}
	i
    }

    pm.readmagicnumber <- function(con) {
	ch <- pm.getc(con)
	if (ch != "P") stop("Not a PNM format file")
	ch <- as.integer(pm.getc(con))
	if (ch < 1 || ch > 6) stop("Unknown PNM format")
	ascii <- FALSE
	if (ch < 4) ascii <- TRUE
	if (ch == 1 || ch == 4) type <- "pbm"
	else if (ch == 2 || ch == 5) type <- "pgm"
	else if (ch == 3 || ch == 6) type <- "ppm"
	res <- list(type=type, ascii=ascii)
	res
    }

    magic <- pm.readmagicnumber(con)

    nc <- pm.getuint(con)
    nr <- pm.getuint(con)
    if (magic$type != "pbm") maxval <- pm.getuint(con)
    else maxval <- 1

    datastart <- seek(con)
    
    seek(con, 0)
    
    if (nc < 0 || nr < 0 || maxval < 1 || maxval > 65535)
	warning(paste("Possible error reading heading: nc:", nc,
	"nr:", nr, "maxval:", maxval))
    
    res <- list(nc = nc, nr = nr, maxval = maxval, type=magic$type,
		datastart=datastart, ascii=magic$ascii)

    invisible(res)
}


