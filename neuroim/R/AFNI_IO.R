#' @import stringr
NULL



#' @keywords internal
parseIntAttribute <- function(line){ 
	as.integer(str_split(str_trim(paste(line, collapse=" ")), "\\s+")[[1]])
}

#' @keywords internal
parseFloatAttribute <- function(line){ 
	as.numeric(str_split(str_trim(paste(line, collapse=" ")), "\\s+")[[1]])
}
	
#' @keywords internal
parseStringAttribute <- function(line) {
	
	res <- str_split(line, "~")[[1]]
	if (length(res) > 1) {
		res[1] <- str_sub(res[1], 2)
		res[1:(length(res)-1)]
	} else {
		str_sub(res[1:(length(res)-1)], 2)	
	}
}
	
#' @keywords internal
parseElement <- function(inputLines) {
	atype <- str_trim(str_split(inputLines[[1]], "=")[[1]])[2]	
	name <- str_trim(str_split(inputLines[[2]], "=")[[1]])[2]	
	count <- str_trim(str_split(inputLines[[3]], "=")[[1]])[2]	
		
	content <- if (atype == "string-attribute") {
		parseStringAttribute(inputLines[4:length(inputLines)])
	} else if (atype == "integer-attribute") {
		parseIntAttribute(inputLines[4:length(inputLines)])
	} else if (atype == "float-attribute") {
		parseFloatAttribute(inputLines[4:length(inputLines)])
	} else {
		stop("unrecognized attribute type ", atype)
	}
	list(type=atype, name=name, count=count, content=content)
	
}

#' readAFNIHeader
#
#' @param fileName the name of the AFNI header file (ending in .HEAD)
#' 
readAFNIHeader <- function(fileName) {
	inputLines <- scan(fileName, what=character(), sep="\n", blank.lines.skip=FALSE)
	idx <- which(unlist(lapply(inputLines, function(lin) lin == ""))) + 1
	lastIdx <- length(inputLines)
	attlen <- diff(c(idx, lastIdx+1))
	attpos <- cbind(idx, attlen)
	
	header <- apply(attpos, 1, function(row) {
		a <- row[[1]]
		b <- row[[2]]
		parseElement(inputLines[a:(a+b-1)])
	})

	names(header) <- unlist(lapply(header, "[[", "name"))
	header
	
}

readNIMLSurfaceHeader <- function(fileName) {
  p <- parse_niml_file(fileName)
  whdat <- which(unlist(lapply(p, "[[", "label")) == "SPARSE_DATA")
  dmat <- do.call(cbind, p[[whdat]]$data)
  whind <- which(unlist(lapply(p, "[[", "label")) == "INDEX_LIST")
  idat <- p[[whind]]$data[[1]]
  list(nodeCount=nrow(dmat), nels=ncol(dmat), 
       label=stripExtension(NIML_SURFACE_DSET, basename(fileName)),
       data=dmat, nodes=idat)
}

parse_niml_element <- function(el) {
  items <- strsplit(el, " ")[[1]]
  if (length(items) > 1) {
    items <- items[items != "" & items != ">"] 
    label <- items[1]
    els <- lapply(items[2:length(items)], function(it) {
      keyval <- strsplit(it, "=")[[1]]
      c(key=keyval[1], val=keyval[2])
    })
    keys <- sapply(els, "[[", "key")
    vals <- lapply(els, "[[", "val")
    names(vals) <- keys
    list(label=label, attr=vals)
  } else if (length(items) == 1) {
    list(label=items[[1]], attr=NULL)
  }
}



read_niml_binary <- function(fconn, meta) {
  dtype <- meta$ni_type
  dtype <- strsplit(as.character(dtype), "\\*")[[1]]
  
  if (length(dtype) == 2) {
    nvols <- as.integer(dtype[1])
    type <- dtype[2]
  } else {
    nvols <- 1
    type <- dtype[1]
  }
  
  type <- switch (type,
                  int="integer",
                  double="double",
                  float="double")
  
  
  nels <- as.integer(meta$ni_dimen)
  
  volseq <- lapply(1:nvols, function(i) {
    readBin(fconn, what=type, size=4, n=nels)
  })
  
  volseq
}

parse_niml_header <- function(fconn) {
  out <- c()
  STATE <- "BEGIN"
  while(TRUE ) {
    ch <- readChar(fconn,1)
    print(ch)
    if (length(ch) == 0) {
      break
    } else if (ch == "<" && STATE == "BEGIN") {
      ## open header
      STATE <- "HEADER"
    } else if (ch == ">" && STATE == "HEADER") {
      STATE <- "END"
      break
    } else {
      out <- c(out, ch)
    }
    
    print(STATE)
  }
  
  out <- paste(out, collapse="")
  out <- gsub("\n", "", out)
  out <- gsub("\"", "", out)
  out <- gsub("/", "", out)
  
  ret <- parse_niml_element(str_trim(out))
  ret
}

parse_niml_next <- function(fconn) {
  header <- parse_niml_header(fconn)
  if (!is.null(header$attr) && (header$label == "SPARSE_DATA" || header$label == "INDEX_LIST")) {
    header$data <- read_niml_binary(fconn, header$attr)
    #while (readChar(fconn,1) != ">") { next }
  }
  out <- c()
  STATE <- "BEGIN"
  lastch <- ""
  while(TRUE ) {
    ch <- readChar(fconn,1)
    print(ch)
    if (length(ch) == 0) {
      break
    } else if (ch == "<" && STATE == "BEGIN") {
      STATE <- "CLOSE_TAG"
    } else if (ch == ">" && STATE == "CLOSE_TAG") {
      ## open header
      STATE <- "END"
      break
    } 
  }
  
  header
  
}

parse_niml_file <- function(fname, maxels=10000) {
  fconn <- file(fname)
  open(fconn, open="rb")
  fsize <- file.info(fname)$size
  out <- list()
  elcount <- 1
  out[[elcount]] <- parse_niml_header(fconn)
  while (seek(fconn, where=NA) < fsize && elcount < maxels) {
    elcount <- elcount + 1
    el <- parse_niml_next(fconn)
    print(el$label)
    out[[elcount]] <- el
    print(seek(fconn, where=NA))
  }
  
  out
}


