## s$Id: read.xls.R 1596 2012-08-22 15:45:22Z warnes $

read.xls <- function(xls, sheet = 1, verbose=FALSE, pattern, 
                     na.strings = c("NA","#DIV/0!"), ...,
                     method=c("csv","tsv","tab"),
                     perl="perl")
{
  con <- tfn <- NULL
  on.exit({ 
	err <- FALSE
	if (inherits(con, "connection")) {
		tryCatch(op <- isOpen(con), error = function(x) err <<- TRUE)
		if (!err && op) close(con)
	}
    if (file.exists(tfn)) file.remove(tfn)
  })

  method <- match.arg(method)
  
  ## expand file path, translating ~ to user's home directory, etc.
  xls <- path.expand(xls)


  ## translate from xls to csv/tsv/tab format (returns name of created file)
  perl <- if (missing(perl))
    findPerl(verbose = verbose)
  else
    findPerl(perl, verbose = verbose)
  
  con <- xls2sep(xls, sheet, verbose=verbose, ..., method=method, perl = perl)
  ## While xls2sep returns a connection, we are better off directly
  ## opening the file, so that R can properly handle the encoding.  So,
  ## just grab the full file path to use later, and close the connection.
  tfn <- summary(con)$description  
  close(con)

  if (missing(pattern))
    {
      if(verbose)
        cat("Reading", method, "file ", dQuote(tfn), "...\n")
      
      if(method=="csv")
        retval <- read.csv(tfn, na.strings=na.strings, ...)
      else if (method %in% c("tsv","tab") )
        retval <- read.delim(tfn, na.strings=na.strings, ...)
      else
        stop("Unknown method", method)
        
      if(verbose)
        cat("Done.\n")
    }
  else {
    if(verbose)
      cat("Searching for lines tfntaining pattern ", pattern, "... ")
    idx <- grep(pattern, readLines(tfn))
    if (length(idx) == 0) {
      warning("pattern not found")
      return(NULL)
    }
   if(verbose)
     cat("Done.\n")
    
    if(verbose)
      cat("Reading", method, "file ", dQuote(tfn), "...\n")

    if(method=="csv")
      retval <- read.csv(tfn, skip = idx[1]-1, na.strings=na.strings, ...)
    else if (method %in% c("tsv","tab") )
      retval <- read.delim(tfn, skip = idx[1]-1, na.strings=na.strings, ...)
    else
      stop("Unknown method", method)

    if(verbose)
      cat("Done.\n")     
  }

  retval
}

