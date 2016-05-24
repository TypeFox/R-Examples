## s$Id: read.xls.R 1423 2010-02-21 17:12:30Z ggrothendieck2 $

installXLSXsupport <- function(perl="perl", verbose=FALSE)
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
  temp.dir <- tempdir()
  ##
  ##

  cmd <- "install_modules.pl"
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
      cat("Attempting to automaticall install Perl libraries to support XLSX (Excel 2007+) file format...\n")
      cat("\n")
    }

  ##

  output <- system(cmd, intern=TRUE)

  ##

  if(verbose) cat("Results: ", output, "\n")

  ##

  if( "XLSX" %in% xlsFormats(perl=perl, verbose=verbose) )
	{
		cat("\nPerl XLSX support libraries successfully installed.\n\n")
		invisible(TRUE)
	}
  else
	{
        	stop("\nUnable to install Perl XLSX support libraries.\n\n")
		invisible(FALSE)
	}
}
