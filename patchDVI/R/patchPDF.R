SweavePDFMiktex <- function( Rnw, main=outputname,  
                             cmd="texify --pdf", 
                             options="--tex-option=-synctex=-1 --tex-option=-interaction=nonstopmode",
                             includedir="--tex-option=-include-directory=",
                             stylepath=FALSE,
                             source.code=NULL,
                             make=1,
                             preview=NULL,
			     patchLog = TRUE,
			     sleep = 0,
                             ...) {
    if (!is.null(source.code) && file.exists(source.code))
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
    cat(consoleLog, sep="\n")
    if(!is.null(status) && status) Sys.sleep(sleep)
    message(patchSynctex(sub("\\.tex$", ".synctex", main, ignore.case = TRUE), patchLog = patchLog))
    if (!is.null(preview)) {
	pdf <- sub("\\.tex$", ".pdf", main, ignore.case = TRUE)
    	cmd <- sprintf(preview, pdf)
    	cat(cmd, "\n")
    	system(cmd, wait=FALSE, invisible=FALSE)
    }
}

SweavePDF <- function( Rnw, main=outputname,
                       texinputs=NULL,
                       source.code=NULL,
                       make=1,
                       links = NULL,
                       preview = NULL,
                       patchLog = TRUE,
                       ... ) {
    if (!is.null(source.code) && file.exists(source.code))
    	try(source(source.code, local=TRUE))
    if (sub(".*\\.tex$", "TeX", Rnw, ignore.case = TRUE) == "TeX") 
    	outputname <- Rnw
    else
    	outputname <- SweaveAll(Rnw, make=make, ...)[1]
    consoleLog <- try(texi2dvi(main, pdf=TRUE, texinputs=texinputs, links=links))
    if (patchLog && !inherits(consoleLog, "try-error")) {
        tempLog <- tempfile(fileext = ".log")
        writeLines(consoleLog, tempLog)
        patchLog(tempLog)
        consoleLog <- readLines(tempLog)
    }
    cat(consoleLog, sep = "\n")
    message(patchSynctex(sub("\\.tex$", ".synctex", main, ignore.case=TRUE), patchLog = patchLog))
    pdf <- sub("\\.tex$", ".pdf", main, ignore.case = TRUE)
    if (!is.null(preview)) {
    	cmd <- sprintf(preview, pdf)
    	cat(cmd, "\n")
    	system(cmd, wait=FALSE, invisible=FALSE)
    }   
}

SweaveDVIPDFM <- function(Rnw, main=outputname,
		       latex = "latex",
		       latexOpts = "-synctex=1 -interaction=nonstopmode",
		       dvipdfm = "dvipdfm",
		       dvipdfmOpts = "",
                       texinputs=NULL,
                       source.code=NULL,
                       make=1,
                       preview = NULL,
                       patchLog = TRUE,
                       ... ) {
    if (!is.null(source.code) && file.exists(source.code))
    	try(source(source.code, local=TRUE))
    if (sub(".*\\.tex$", "TeX", Rnw, ignore.case = TRUE) == "TeX") 
    	outputname <- Rnw
    else
    	outputname <- SweaveAll(Rnw, make=make, ...)[1]
    
    otexinputs <- Sys.getenv("TEXINPUTS", unset = NA)
    if(is.na(otexinputs)) {
        on.exit(Sys.unsetenv("TEXINPUTS"))
	    otexinputs <- "."
    } else on.exit(Sys.setenv(TEXINPUTS = otexinputs))
    Sys.setenv(TEXINPUTS = paste(otexinputs, Rtexinputs(), "", sep = .Platform$path.sep))
    cmd <- paste(shQuote(latex), latexOpts, shQuote(main))
    consoleLog <- try(system(cmd, intern = TRUE))
    if (patchLog && !inherits(consoleLog, "try-error")) {
        tempLog <- tempfile(fileext = ".log")
        writeLines(consoleLog, tempLog)
        patchLog(tempLog)
        consoleLog <- readLines(tempLog)
    }
    cat(consoleLog, sep = "\n")
    dvi <- sub("\\.tex$", ".dvi", main, ignore.case = TRUE)
    message(patchSynctex(sub("\\.tex$", ".synctex", main, ignore.case=TRUE), fromDVI = dvi, 
            patchLog = patchLog))
    cmd <- paste(shQuote(dvipdfm), dvipdfmOpts, shQuote(dvi))
    system(cmd)
    if (!is.null(preview)) {
        pdf <- sub("\\.tex$", ".pdf", main, ignore.case = TRUE)
    	cmd <- sprintf(preview, pdf)
    	cat(cmd, "\n")
    	system(cmd, wait=FALSE, invisible=FALSE)
    }   
}

rawToLines <- function(raw) {
    temp <- tempfile()
    on.exit(unlink(temp))
    writeBin(raw, temp)
    readLines(temp, warn=FALSE)
}

pdfEOF <- function(con) {
    seek(con, -1024, "end")
    tail <- rawToChar(readBin(con, "raw", 1024), multiple=TRUE)
    stop <- rev(which(tail == "F"))
    for (last in stop) {
	if (all(tail[last - 4:0] == c("%", "%", "E", "O", "F")))
	    return( last - 4 + seek(con, 0) - 1024 - 1 )
    }	    
    stop("%%EOF marker not found")
}

pdfStartxrefs <- function(con, eof=pdfEOF(con)) {
    seek(con, eof - 20)
    tail <- rawToLines(readBin(con, "raw", 20))
    startxref <- which(tail == "startxref")
    as.numeric(tail[startxref+1])
}

pdfXrefsCompressed <- function(pdfname) {
    con <- file(pdfname, "rb")
    on.exit(close(con))
    seek(con, pdfStartxrefs(con, pdfEOF(con)))
    xrefs <- rawToLines(readBin(con, "raw", 50))
    grepl(" obj <<$", xrefs[1])
}

pdfXrefblock <- function(con, start=pdfStartxrefs(con, eof), eof=pdfEOF(con)) {
    seek(con, start)
    xrefs <- rawToLines(readBin(con, "raw", eof - start))
    trailer <- which(xrefs == "trailer")[1]
    if (is.na(trailer)) stop("PDF file can't be read--may be using compression.")
    tail <- xrefs[trailer:length(xrefs)]
    xrefs <- xrefs[1:(trailer-1)]
    line <- 2
    offsets <- numeric(0)
    generations <- numeric(0)
    free <- logical(0)
    maxobj <- 0
    while (line <= length(xrefs)) {
	  head <- scan(temp <- textConnection(xrefs[line]), 
	               what=list(first=0, count=0),
	               quiet = TRUE)
	  close(temp)
	  newmax <- head$first + head$count 
	  if (newmax > maxobj) {
	      length(offsets) <- newmax
	      length(generations) <- newmax
	      length(free) <- newmax
	  }
	  if (head$count > 0) {
	      body <- scan(temp <- textConnection(xrefs[line + 1:head$count]),
	                   what = list(offsets=0, generations=0, free="n"),
	                   quiet = TRUE)
	      close(temp)
	      inds <- head$first + 1:head$count
	      offsets[inds] <- body$offsets
	      generations[inds] <- body$generations
	      free[inds] <- body$free == "f"
	  }
	  line <- line + head$count + 1
    }
    # Now clean up the tail dictionary
    stop <- which(tail == "startxref")[1]
    tail <- tail[2:(stop-1)]
    tail[1] <- sub("^<< *", "", tail[1])
    tail[length(tail)] <- sub(" *>>$", "", tail[length(tail)])    
    list(xrefs = data.frame(offsets = offsets, generations = generations, free=free),
         dict = tail)
         
}

dictVal <- function(dict, key) {
    dict <- unlist(strsplit(paste(dict, collapse=" "), "/"))
    dict <- grep(paste("^", key, sep=""), dict, value=TRUE)
    return(substr(dict, nchar(key)+1, nchar(dict)))
}

pdfXrefblocks <- function(con, collapse = TRUE) {
    result <- list()
    eof <- pdfEOF(con)
    start <- pdfStartxrefs(con, eof)
    rows <- NULL
    repeat {
        block <- pdfXrefblock(con, start, eof)
        if (collapse) {
            if (is.null(rows)) {
		rows <- as.numeric(dictVal(block$dict, "Size"))
		NAs <- rep(NA_real_, rows)
		xref <- data.frame(offsets = NAs, generations = NAs, free= as.logical(NAs))
	    }
            replace <- ((1:rows) <= nrow(block$xrefs)) & is.na(xref$offsets)
            xref[replace,] <- block$xrefs[replace,]
       	} else
            result <- c(result, list(block))
        start <- as.numeric(dictVal(block$dict, "Prev"))
        if (!length(start)) break
    }
    if (collapse) 
        return(xref)
    else
        return(result)
}

pdfFindobj <- function(con, pattern) {
    xrefs <- pdfXrefblocks(con)
    xrefs <- subset(xrefs, !xrefs$free)
    o <- order(xrefs$offsets)
    xrefs <- xrefs[o,]
    result <- character(0)
    for (i in 1:nrow(xrefs)) {
    	seek(con, xrefs$offsets[i])
    	head <- readChar(con, nchar(pattern) + 30)
    	head <- sub("^[[:space:][:digit:]]*obj[[:space:]]*", "", head)
    	if (grepl(pattern, head)) {
            seek(con, xrefs$offsets[i])
            obj <- rawToLines(readBin(con, "raw", xrefs$offsets[i+1]-xrefs$offsets[i]))
            # FIXME:  this will match any endobj, not just the one we want
            stop <- grep("endobj", obj)
            obj <- paste(obj[1:stop], collapse=" ")
            obj <- sub("^[[:space:][:digit:]]*obj[[:space:]]*", "", obj)
            obj <- sub("[[:space:]]*endobj.*", "", obj)
            result <- c(result, obj)
        }
    }
    result
}

pdfobjs <- function(file, pattern) {
    if (is.character(file)) {
        file <- file(file, "rb")
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    if (!isOpen(file, "rb")) {
        open(file, "rb")
        on.exit(close(file))
    }
    pdfFindobj(file, pattern)
}

pdfStreams <- function(file, pattern) {
    streamhead <- 
    streams <- pdfobjs(file, "<<\n/Length[[:space:]]+[[:digit:]]+[[:space:]]*\n>>\nstream\n")
    streams <- sub("^<<[[:space:]]*/Length[[:space:]]+[[:digit:]]+[[:space:]]*>>[[:space:]]*stream[[:space:]]*", "", streams)
    streams <- sub("[[:space:]]*endstream$", "", streams)
    streams <- grep(pattern, streams, value=TRUE)
    streams
}

syncFiles <- function(lines) {
    inputs <- grep("^Input:", lines)
    filenames <- sub("^Input:[[:digit:]]*:","",lines[inputs])
    filepaths <- dirname(filenames)
    filenames <- basename(filenames)
    nodot <- grep("^[^.]*$", filenames)
    filenames[nodot] <- paste(filenames[nodot], ".tex", sep="")
    filenums <- sub("^Input:","",lines[inputs])
    filenums <- as.numeric(sub(":.*","",filenums))
    o <- order(filenums)
    data.frame(tag=filenums[o], path=filepaths[o], name=filenames[o])
}

parseConcords <- function(lines) {
   parseConcord <- function(split) {
    	oldname <- split[2]
    	newsrc <- split[3]
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
    	list(oldname=oldname, newsrc=newsrc, concord=concord, ofs=ofs)
    }
    concords <- strsplit(lines, ":")
    concords <- lapply(concords, parseConcord)
    names(concords) <- sapply(concords, function(x) x$oldname)
    concords
}

grepConcords <- function(pdfname) {
    # grepRaw is new in R 2.13.0
    if (!exists("grepRaw")) return(character(0))
    
    size <- file.info(pdfname)$size
    if (is.na(size)) stop(pdfname, " not found")

    buffer <- readBin(pdfname, "raw", size)
    result <- grepRaw("concordance:[^:\n[:space:]]+:[^:\n[:space:]]+:(ofs [[:digit:]]+:)?[[:digit:]][-[:digit:] ]*", 
            buffer, fixed=FALSE, all=TRUE, value=TRUE)
    if (!length(result)) character(0)
    else sapply(result, rawToChar)
}
    
patchSynctex <- function(f, newname=f, uncompress="pdftk %s output %s uncompress",
		         fromDVI = NULL, patchLog = TRUE) {
    basename <- tools::file_path_sans_ext(f)
    compressed <- FALSE
    if (!file.exists(f)) {
    	f <- paste(f, ".gz", sep="")
    	if (file.exists(f)) {
    	    compressed <- TRUE
    	    force(newname)
    	    f <- gzfile(f)
    	}
    }
    lines <- try(readLines(f, warn=FALSE), silent=TRUE)
    if (inherits(lines, "try-error")) 
    	return(paste(f,"cannot be read, no patching done."))

    files <- syncFiles(lines)
    if (is.null(fromDVI)) {
	pdfname <- file.path(files$path[1], paste(sub(".tex", "", files$name[1]), ".pdf", sep=""))
    
	concords <- parseConcords(grepConcords(pdfname))
	if (!length(concords))
	   concords <- readConcords(list.files(as.character(files$path[1]),
	   				       pattern = "-concordance.tex$"))
	
	if (!length(concords)) {
	    if (pdfXrefsCompressed(pdfname)) {
		if (missing(uncompress) && nchar(Sys.which("pdftk")) == 0) 
		    return(paste("No patches made:\n", pdfname, 
    	                 "\nis compressed and no uncompressor was found."))
		oldname <- pdfname
		pdfname <- tempfile()
		system(sprintf(uncompress, oldname, pdfname))
	    }
    
	    concords <- parseConcords(pdfStreams(pdfname, "^concordance:"))
	    if (!length(concords)) # try older buggy format
		concords <- parseConcords(pdfobjs(pdfname, "^concordance:"))
	}
    } else {
	specials <- DVIspecials(fromDVI)
	concordind <- grep("^concordance:", specials)
	concords <- parseConcords(specials[concordind])
	if (length(concords)) {
	    specials[concordind] <- NA
	  #  setDVIspecials(fromDVI, specials)
	}
    }
    if (patchLog)
    	patchLog(paste0(basename, ".log"), concords = concords)
    re <- "^([vhxkgr$[(])([[:digit:]]+),([[:digit:]]+)([^[:digit:]].*)"
    srcrefind <- grep(re, lines)
    srcrefs <- lines[srcrefind]
    
    ops <- sub(re, "\\1", srcrefs)
    tags <- sub(re, "\\2", srcrefs)
    linenums <- sub(re, "\\3", srcrefs)
    rest <- sub(re, "\\4", srcrefs)
    
    linenums <- as.integer(linenums)
     
    changed <- rep(FALSE, length(tags))
    newtags <- c()
    maxtag <- max(files$tag)
    for (i in seq_along(concords)) {
        n <- names(concords)[i]
        ofs <- concords[[i]]$ofs
        concord <- concords[[i]]$concord
        newsrc <- concords[[i]]$newsrc
        if (!(newsrc %in% names(newtags))) {
            maxtag <- maxtag + 1
            newtags <- c(newtags, maxtag)
            names(newtags)[length(newtags)] <- newsrc
        }
        tag <- files$tag[files$name == n]
        if (length(tag) == 1) {
    	    subset <- (tags == tag) & (linenums > ofs) & (linenums <= ofs + length(concord))
    	    linenums[subset] <- concord[linenums[subset] - ofs]
    	    tags[subset] <- newtags[newsrc]
    	    changed[subset] <- TRUE
    	}
    }
    
    if (any(changed)) {
    	newrefs <- ifelse(changed, paste(ops, tags, ",", linenums, rest, sep=""), srcrefs)
    	lines[srcrefind] <- newrefs
    	firstInput <- grep("^Input:", lines)[1]
    	lines <- c(lines[1:firstInput],
    	           paste("Input:", newtags, ":", names(newtags), sep=""),
    	           lines[(firstInput+1):length(lines)])
    }
    # Texworks on Windows assumes Unix line endings; make sure
    # we create those.
    con <- if (compressed) gzfile(newname, "wb") else file(newname, "wb")
    on.exit(close(con))
    writeLines(lines, con, sep="\xa")
    changes <- sum(changed) + length(newtags)
    msg <- paste(changes, "patches made.") 
    if (!changes)
    	msg <- paste(msg, "Did you set \\SweaveOpts{concordance=TRUE}?")
    msg
}