SweaveMiktex <- function(Rnw, 
                         main=outputname, 
                         cmd="texify",
                         options="--tex-option=-src-specials --tex-option=-interaction=nonstopmode",
                         includedir="--tex-option=--include-directory=",
                         stylepath=FALSE,
                         source.code=NULL,
                         make=1,
                         preview='yap "%s"',
			 patchLog = TRUE,
			 sleep = 0,
                         ...) {
    if (!is.null(source.code))
    	try(source(source.code, local=TRUE))
    if (sub(".*\\.tex$", "TeX", Rnw, ignore.case = TRUE) == "TeX") 
    	outputname <- Rnw
    else
    	outputname <- SweaveAll(Rnw, make=make, stylepath=stylepath, ...)[1]
    cmd <- paste(cmd, " ", options, " ", includedir, Rtexinputs(),
                 " ", main, sep="")    	
    cat(cmd, "\n")
    consoleLog <- try(system(cmd, intern = TRUE))
    status <- attr(consoleLog, "status")
    if (patchLog && !inherits(consoleLog, "try-error")) {
        tempLog <- tempfile(fileext = ".log")
        writeLines(consoleLog, tempLog)
        patchLog(tempLog)
        consoleLog <- readLines(tempLog)
    }
    cat(consoleLog, sep = "\n")
    if (!is.null(status) && status) Sys.sleep(sleep)
    dvi <- sub("\\.tex$", ".dvi", main, ignore.case = TRUE)
    message(patchDVI(dvi, patchLog = patchLog))
    if (!is.null(preview)) {
    	cmd <- sprintf(preview, dvi)
    	cat(cmd, "\n")
    	system(cmd, wait=FALSE, invisible=FALSE)
    }
}

SweaveDVI <- function( Rnw, main=outputname,
                       texinputs=NULL,
                       source.code=NULL,
                       make=1,
                       links=NULL,
                       preview=NULL,
		       patchLog = TRUE,
                       ... ) {
    if (!is.null(source.code) && file.exists(source.code))
    	try(source(source.code, local=TRUE))
    if (sub(".*\\.tex$", "TeX", Rnw, ignore.case = TRUE) == "TeX") 
    	outputname <- Rnw
    else
    	outputname <- SweaveAll(Rnw, make=make, ...)[1]
    consoleLog <- try(texi2dvi(main, pdf=FALSE, texinputs=texinputs, links=links))
    if (patchLog && !inherits(consoleLog, "try-error")) {
        tempLog <- tempfile(fileext = ".log")
        writeLines(consoleLog, tempLog)
        patchLog(tempLog)
        consoleLog <- readLines(tempLog)
    }
    cat(consoleLog, sep = "\n")
    dvi <- sub("\\.tex$", ".dvi", main, ignore.case=TRUE)
    message(patchDVI(dvi, patchLog = patchLog))
    if (!is.null(preview)) {
    	cmd <- sprintf(preview, dvi)
    	cat(cmd, "\n")
    	system(cmd, wait=FALSE, invisible=FALSE)
    }    
}

readDVI <- function(f, show=c("bop", "special", "fntdef", "preamble")) {
    size <- file.info(f)$size
    con <- file(f, "rb")
    bytes <- readBin(con, "raw", size)
    close(con)
    
    parmsizes <- c( rep(0, 128), 1:4, 8, 1:4, 8,   # 0 to 137
                    0, 44, 0,0,0, 1:4, 0, 1:4, 0,  # 138 to 152
                    1:4, 1:4, 0, 1:4, 0, 1:4, 	   # 153 to 170
                    rep(0, 64), 1:4, 1:4, 15:18,   # 171 to 246
                    14, 28, 5, 9, rep(NA, 5))      # 247 to 255
                    
    pos <- 0
    opcode <- 0
    while (pos < size) {
        pos <- pos+1
    	opcode <- as.integer(bytes[pos])
    	parmsize <- parmsizes[opcode + 1]
    	
    	if (opcode < 139L) { # do nothing
    	} else if (opcode == 139L && "bop" %in% show) {  # bop
    	    counters <- readBin(bytes[pos + 1:(parmsize-4)], "integer", n= 10, size=4, endian="big")
    	    prev <- readBin(bytes[pos + 40 + 1:4], "integer", size=4, endian="big")
    	    cat("bop at pos=", pos, ": counters=", paste(counters, collapse=" "), " prev=", prev, "\n")
    	} else if (opcode < 239L) { # do nothing 
    	} else if (opcode < 243L) {   # xxxi
    	    k <- readBin(bytes[pos + (1:parmsize)], "integer", size=parmsize, endian="big")
    	    if ("special" %in% show) {
    	    	special <- readChar(bytes[(pos + parmsize + 1):(pos + parmsize + k)], k)
    	    	cat("special at pos=", pos," is ", special, "\n")
    	    }
    	    parmsize <- parmsize + k
    	} else if (opcode < 247L) { # fnt def i
    	    a <- as.integer(bytes[pos + parmsize ])
    	    if ("fntdef" %in% show) {
     	    	fntname <- readChar(bytes[(pos + parmsize + 1):(pos + parmsize + a)], a)
   	    	cat("fnt def at pos=", pos, " is ", fntname, "\n")
   	    }
    	    parmsize <- parmsize + a 
    	} else if (opcode == 247L) {       # pre
    	    k <- as.integer(bytes[pos + parmsize])
	    if ("preamble" %in% show) {
    	    	comment <- readChar(bytes[(pos + parmsize + 1):(pos + parmsize + k)], k)        	    
    	    	cat("preamble at pos=", pos, " with comment ", comment, "\n")
	    }
    	    parmsize <- parmsize + k
    	} else if (opcode == 249L) break
	pos <- pos + parmsize
    }
}

DVIspecials <- function(f) {
    size <- file.info(f)$size
    con <- file(f, "rb")
    bytes <- readBin(con, "raw", size)
    close(con)
    .Call(dviSpecials, bytes)
}

setDVIspecials <- function(f, newspecials, newname=f) {
    size <- file.info(f)$size
    con <- file(f, "r+b")
    on.exit(close(con))
    bytes <- readBin(con, "raw", size)
    bytes <- .Call(setDviSpecials, bytes, as.character(newspecials))
    close(con)
    con <- file(newname, "wb")
    writeBin(bytes, con)
}

patchDVI <- function(f, newname=f, patchLog = TRUE) {
    specials <- DVIspecials(f)
    
    concordind <- grep("^concordance:", specials)
    concords <- specials[concordind]
    
    parseConcord <- function(split) {
    	oldname <- split[2]
    	newname <- split[3]
    	if (length(split) == 4) {
    	    ofs <- 0
    	    vi <- 4
    	} else {
    	    ofs <- as.integer(sub("^ofs ([0-9]+)", "\\1", split[4]))
    	    vi <- 5
    	}
    	values <- as.integer(strsplit(split[vi], " ")[[1]])
    	firstline <- values[1]
    	rledata <- matrix(values[-1], nrow=2)
    	rle <- structure(list(lengths=rledata[1,], values=rledata[2,]), class="rle")
    	diffs <- inverse.rle(rle)
	concord <- c(firstline, firstline + cumsum(diffs))
    	list(oldname=oldname, newname=newname, concord=concord, ofs=ofs)
    }
    concords <- strsplit(concords, ":")
    concords <- lapply(concords, parseConcord)
    names(concords) <- sapply(concords, function(x) x$oldname)

    if (patchLog)
    	patchLog(paste0(tools::file_path_sans_ext(f), ".log"), concords = concords)

    srcrefind <- grep("^src:", specials)
    srcrefs <- specials[srcrefind]
    
    linenums <- sub("^src:([0-9]+).*$", "\\1", srcrefs)
    filenames <- substr(srcrefs, 5+nchar(linenums), 10000)
    # \input{file} doesn't put .tex on it
    noext <- !grepl("\\.", filenames)
    filenames[noext] <- paste(filenames[noext], ".tex", sep="")
    filenames <- myNormalizePath(filenames)
    linenums <- as.integer(linenums)
    
    changed <- rep(FALSE, length(filenames))
    for (i in seq_along(concords)) {
    	n <- names(concords)[i]
    	ofs <- concords[[i]]$ofs
    	concord <- concords[[i]]$concord
    	subset <- (filenames == myNormalizePath(n)) & (linenums > ofs) & (linenums <= ofs + length(concord))
    	linenums[subset] <- concord[linenums[subset] - ofs]
    	filenames[subset] <- concords[[i]]$newname
    	changed[subset] <- TRUE
    }
    
    newrefs <- ifelse(changed, paste("src:", linenums, filenames, sep=""), srcrefs)
    
    specials[srcrefind] <- newrefs
    specials[concordind] <- NA
    
    if (any(changed) || length(concordind))
    	setDVIspecials(f, specials, newname)
    
    changes <- sum(changed) + length(concordind)
    msg <- paste(changes, "patches made.")
    if (!changes)
    	msg <- paste(msg, "Did you set \\SweaveOpts{concordance=TRUE}?")
    msg
}

