## s$Id: read.xls.R 1423 2010-02-21 17:12:30Z ggrothendieck2 $

xlsFormats <- function(perl="perl", verbose=FALSE)
{
  ## determine proper path to perl executable
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

  cmd <- "supportedFormats.pl"
  sc <- file.path(perl.dir, cmd)
  
  ##
  ##

  ##
  ## execution command

  cmd <- paste(shQuote(perl), shQuote(sc), sep=" ")

  ##

  if(verbose)
    {
      cat("\n")
      cat("Determining supported formats...\n")
      cat("\n")
    }

  ##

  output <- system(cmd, intern=TRUE)

  ##

  if(verbose) cat("Results: ", output, "\n")

  ##

  retval <- unlist( strsplit(output," "))
  retval <- retval[ -c(1,2) ]

  return(retval)
}
