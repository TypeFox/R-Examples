## s$Id: xls2sep.R 1639 2013-01-14 20:47:57Z warnes $

xls2csv <- function(xls, sheet=1, verbose=FALSE, blank.lines.skip=TRUE,
                    ..., perl="perl")
  xls2sep(xls=xls, sheet=sheet, verbose=verbose,
          blank.lines.skip=blank.lines.skip, ..., method="csv",
          perl=perl)

xls2tab <- function(xls, sheet=1, verbose=FALSE, blank.lines.skip=TRUE,
                    ..., perl="perl")
  xls2sep(xls=xls, sheet=sheet, verbose=verbose,
          blank.lines.skip=blank.lines.skip, ..., method="tab",
          perl=perl) 

xls2tsv <- function(xls, sheet=1, verbose=FALSE, blank.lines.skip=TRUE,
                    ..., perl="perl")
  xls2sep(xls=xls, sheet=sheet, verbose=verbose,
          blank.lines.skip=blank.lines.skip, ..., method="tsv",
          perl=perl) 

xls2sep <- function(xls,
                    sheet=1,
                    verbose=FALSE,
                    blank.lines.skip=TRUE,
                    ...,
                    method=c("csv","tsv","tab"),
                    perl = perl)
  {
    
    method <- match.arg(method)
    
	perl <- if (missing(perl))
		findPerl(verbose = verbose)
	else
		findPerl(perl, verbose = verbose)

    ##
    ## directories
    package.dir <- find.package('gdata')
    perl.dir <- file.path(package.dir,'perl')
    ##

    ##
    ## filesheet
    tf <- NULL
    if ( substring(xls, 1, 7) == "http://" ||
         substring(xls, 1, 6) == "ftp://" )
      {
        tf <- paste(tempfile(), "xls", sep = ".")
        if(verbose)
          cat("Downloading",
              dQuote(xls), " to ",
              dQuote(tf), "...\n")
        download.file(xls, tf, mode = "wb")
        if(verbose) cat("Done.\n")
        xls <- tf
      }

    if(method=="csv")
      {
        script <- file.path(perl.dir,'xls2csv.pl')
        targetFile <- paste(tempfile(), "csv", sep = ".")
      }
    else if(method=="tab")
      {
        script <- file.path(perl.dir,'xls2tab.pl')
        targetFile <- paste(tempfile(), "tab", sep = ".")
      }
    else if(method=="tsv")
      {
        script <- file.path(perl.dir,'xls2tsv.pl')
        targetFile <- paste(tempfile(), "tsv", sep = ".")
      }
    else
      {
        stop("Unknown method", method)
      }
    
    ##
    ##

    ##
    ## blank.lines.skip
    ##
    if (blank.lines.skip)
      skipBlank=""
    else
      skipBlank="-s"

    ##
    ## execution command
    cmd <- paste(shQuote(perl),
                 shQuote(script),
                 skipBlank,  # flag is not quoted
                 shQuote(xls),
                 shQuote(targetFile),
                 shQuote(sheet),
                 sep=" ")
    ##
    ##

    if(verbose)
      {
        cat("\n")
        cat("Converting xls file\n")
        cat("   ", dQuote(xls), "\n")
        cat("to", method, " file \n")
        cat("   ", dQuote(targetFile), "\n")
        cat("... \n\n")
      }
    
    ##
    ## do the translation
    if(verbose)  cat("Executing '", cmd, "'... \n\n")

    results <- try(system(cmd, intern=!verbose))

    if(inherits(results, "try-error"))
      stop( "Unable to read xls file '", xls, "':", results )
    
    if(verbose) cat(results,"\n\n")
    if (verbose) cat("Done.\n\n")

    ##
    ## check that the target file was created
    ##
    if(!file.exists(targetFile))
      stop( "Intermediate file '", targetFile, "' missing!" )

    ## Creae a file object to hand to the next stage..
    retval <- try(file(targetFile))
    if(inherits(retval, "try-error"))

      stop("Unable to open intermediate file '", targetFile, "':",
           retval)
    
    return(retval)
    
  }

