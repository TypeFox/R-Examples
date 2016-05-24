dataURI <- function(data, mime="", encoding="base64", file) {
   if (!is.null(encoding) && !isTRUE(encoding == "base64")) stop('encoding must be either NULL or "base64"') 
   prefix <- paste("data:", as.character(mime)[1], if (!is.null(encoding)) ";base64", ",", sep ='')
   if (!missing(file)) {
      if (!missing(data)) stop("data and file are mutually exclusive")
      data <- con <- file(file, "rb")
      on.exit(close(con))
   }
   if (inherits(data, "connection")) {
      if (isTRUE(summary(data)$text == "binary")) {
         l <- list()
	 while (length(r <- readBin(data, raw(0), 1048576L))) l <- c(l, r)
	 data <- unlist(l)
      } else data <- readLines(data)
   }
   if (!is.raw(data)) data <- paste(as.character(data), collapse='\n')
   paste(prefix, if (is.null(encoding)) .Call(C_URIencode, data, NULL) else .Call(B64_encode, data, 0L, NULL), sep='')
}
