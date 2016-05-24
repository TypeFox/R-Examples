env2Cdf <- function(env, celfile, overwrite=FALSE, verbose=TRUE, ...) {
  require(env, character.only=TRUE) || stop("Package not loaded: ", env)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  env2List <- function(u) {
    pmi <- u[,"pm"];
    mmi <- u[,"mm"];

    x <- (u-1) %% nrows;     # <= BUG?!? /HB 2010-05-20
    y <- (u-1) %/% nrows;    # <= BUG?!? /HB 2010-05-20
    nr <- nrow(u);

    if (all(is.na(mmi))) {
      # do PM only
      # so far this is not implemented
      #o <- which(!is.na(x));
      #g <- list(list(x=u$x[o], y=u$y[o],
      #       pbase=rep("T", times=nr), tbase=rep("A", times=nr),
      #       atom=v, indexpos=v, groupdirection="sense",
      #       natoms=nr, ncellsperatom=1));
      #names(g) <- id;
      #list(unittype=1, unitdirection=1, groups=g, natoms=nr, ncells=nr,
      #   ncellsperatom=1, unitnumber=1);
      return(NULL);
    } else {
      # do PM and MM
      v <- 0:(nr-1);
      atom <- rep(v, times=2);
      o <- order(atom);
      g <- list(list(x=c(x)[o], y=c(y)[o],
             pbase=rep(c("A","T"), times=nr),
             tbase=rep(c("T","T"), times=nr),
             atom=atom[o], indexpos=atom[o], groupdirection="sense",
             natoms=nr, ncellsperatom=2)
           );
      #names(g) <- id;
      list(unittype=1, unitdirection=1, groups=g, natoms=nr, ncells=nr*2,
           ncellsperatom=2, unitnumber=1);
    } # if (...)
  } # env2List()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  if (verbose) {
    cat("Reading environment: ", env, ".\n", sep="");
  }
  ffs <- as.list(get(env));

  if (verbose) {
    cat("Reading CEL file header.\n");
  }

  celHead <- .readCelHeader(celfile);
  nrows <- celHead$rows;
  ncols <- celHead$rows;  # <= BUG?!? /HB 2010-05-20

  nunits <- length(ffs);

  ### make the input-list for the writeCdf-function
  if (verbose) {
    cat("Creating CDF list for ", nunits, " units.\n", sep="");
  }

  newCdfList <- lapply(ffs, FUN=env2List);
  #return(newCdfList)

  if (verbose) {
    cat("Adding group names to CDF list.\n");
  }

  nm <- names(newCdfList);
  for (ii in 1:nunits) {
    names( newCdfList[[ii]]$groups ) <- nm[ii];
    newCdfList[[ii]]$unitnumber <- ii;
  }

  pdName <- gsub("cdf", "", env, fixed=TRUE);

  ## creating the cdf-header;
  pathname <- sprintf("%s.cdf", pdName);
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nunits, nqcunits=0,
    refseq="", chiptype=pdName, filename=pathname, rows=nrows, cols=ncols,
    probesets= nunits, qcprobesets=0, reference="");

  ### writing the cdf-file (binary-file)
  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  res <- .writeCdf(pathnameT, cdfheader=newCdfHeader, cdf=newCdfList,
           cdfqc=NULL, overwrite=overwrite, verbose=verbose);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  invisible(pathname);
} # env2Cdf()

# For backward compatibility
Env2Cdf <- env2Cdf


############################################################################
# HISTORY:
# 2011-08-30 [HB]
# o Now env2Cdf() returns the pathname to the written CDF.
# 2010-09-29 [HB]
# o ROBUSTNESS: Now the writing of the CDF file is atomic by first writing
#   to a temporary file which is then renamed.
# 2010-05-20 [HB]
# o Renamed Env2Cdf() to env2Cdf().  Keeping old one for backward
#   compatibility.
# o Some basic code clean up.
# 2009-01-13 [MR]
# o Added. "This script has been written to generate a .cdf-file from an
#   "XXXXcdf" package, such as the Bioconductor 'metadata' packages.
#   The original was written by Samuel Wuest, modified by Mark Robinson
#   (around 12 Jan 2009) to be generic."
# 2008-??-?? [SW]
# o Created by Samuel Wuest.
############################################################################
