###############################################################################
##
## testPerl.R
##
## Test for the presence of Perl and the required modules
##
## Copyright 2015, Marc Schwartz <marc_schwartz@me.com>
##
## This software is distributed under the terms of the GNU General
## Public License Version 2, June 1991.  

testPerl <- function(perl = "perl", verbose = TRUE) {
  require(WriteXLS)

  ## Get path to WriteXLS Perl tree
  Perl.Path <- system.file("Perl", package = "WriteXLS")

  ## Check For Perl first
  res <- Sys.which(perl)
  
  if ((res == "") & (verbose)) {
    message("\nPerl was not found on your system. Either check $PATH if installed or please install Perl.\n",
            paste("For more information see: ", system.file('INSTALL', package='WriteXLS')), "\n")
    
    invisible(FALSE)
  } else {
    if (verbose)
      message("Perl found.\n")

    PerlModules <- c("Archive/Zip.pm",
                     "OLE/Storage_Lite.pm",
                     "Parse/RecDescent.pm",
                     "Getopt/Long.pm",
                     "File/Basename.pm",
                     "Spreadsheet/WriteExcel.pm",
                     "Excel/Writer/XLSX.pm",
                     "Text/CSV_PP.pm")

    Found <- rep(FALSE, length(PerlModules))
  
    ## Get Perl's current @INC array to search for required modules
    CMD <- paste(perl, "-e \"print join('\n', @INC);\"")
    PERLINC <- system(CMD, intern = TRUE)

    ## Add WriteXLS Perl tree
    ## Substitute for current directory "." in PERLINC if there
    ## Takes a long time to search and modules unlikely to be there anyway
    CurrDir <- grep("^\\.$", PERLINC)
    if (length(CurrDir) > 0) {
      PERLINC[CurrDir] <- Perl.Path
    } else {
      PERLINC <- c(PERLINC, Perl.Path)
    }

    ## Get rid of non-existing or non-readable paths
    PERLINC.exists <- PERLINC[file_test("-d", PERLINC)]

    ## glob the remaining paths and perl modules
    Level1 <- Sys.glob(file.path(PERLINC.exists, "*.pm"))
    Level2 <- Sys.glob(file.path(PERLINC.exists, "*", "*.pm"))
    Level3 <- Sys.glob(file.path(PERLINC.exists, "*", "*", "*.pm"))    
    L1L2L3 <- c(Level1, Level2, Level3)

    ## Search for perl modules
    for (i in seq(along = PerlModules)) {
      Found[i] <- any(grep(PerlModules[i], L1L2L3))
    }

    if (!all(Found)) {
      if (verbose) {
        Missing <- paste(sub("\\.pm", "", gsub("/", "::", PerlModules[!Found])), collapse = "\n")
        message("The following Perl modules were not found on this system:\n")
        message(Missing, "\n")
        message("If you have more than one Perl installation, be sure the correct one was used here.\n")
        message("Otherwise, please install the missing modules. ", paste("For more information see: ", system.file('INSTALL', package='WriteXLS')), "\n")
      }
      invisible(FALSE)
    } else {
      if (verbose)
        message("All required Perl modules were found.\n")
      invisible(TRUE)
    }
  }
}
