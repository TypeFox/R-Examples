if (!require("WriteXLS")) stop("Package WriteXLS is not available.")

z <- Sys.getenv("_R_CHECK_HAVE_PERLCSVXS_")

Sys.info()

if(identical(as.logical(z), TRUE)) { 
  if (!testPerl(verbose=FALSE)) 
     stop("Perl or required modules used by WriteXLS are not functional.") 
  } else {
     message("PERL or modules not available. Skipping tests.\n")
     message("_R_CHECK_HAVE_PERLCSVXS_ setting ", z, "\n")
     }
