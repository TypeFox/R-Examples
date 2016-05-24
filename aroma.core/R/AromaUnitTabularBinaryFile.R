###########################################################################/**
# @RdocClass AromaUnitTabularBinaryFile
#
# @title "The AromaUnitTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  A AromaUnitTabularBinaryFile is an @see "AromaTabularBinaryFile" with
#  the constraint that the rows map one-to-one to, and in the same order as,
#  the units in a annotation chip type file (e.g. CDF file).  
#  The (full) chip type of the annotation chip type file is given by the
#  mandatory file footer \code{chipType}.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# %\seealso{
# % @see "AromaCellTabularBinaryFile".
# %}
#*/########################################################################### 
setConstructorS3("AromaUnitTabularBinaryFile", function(...) {
  extend(AromaMicroarrayTabularBinaryFile(...), c("AromaUnitTabularBinaryFile", uses("UnitAnnotationDataFile")),
    "cached:.unf" = NULL
  );
})


setMethodS3("nbrOfUnits", "AromaUnitTabularBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("byChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, nbrOfUnits=NULL, validate=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'nbrOfUnits':
  if (!is.null(nbrOfUnits)) {
    nbrOfUnits <- Arguments$getInteger(nbrOfUnits, range=c(0,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for all possible matches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnames <- findByChipType(static, chipType=chipType, tags=tags, 
                                                     firstOnly=FALSE, ...);
  if (is.null(pathnames)) {
    ext <- getDefaultExtension(static);
    note <- attr(ext, "note");
    msg <- sprintf("Failed to create %s object. Could not locate an annotation data file for chip type '%s'", class(static)[1], chipType);
    if (is.null(tags)) {
      msg <- sprintf("%s (without requiring any tags)", msg);
    } else {
      msg <- sprintf("%s with tags '%s'", msg, paste(tags, collapse=","));
    }
    msg <- sprintf("%s and with filename extension '%s'", msg, ext);
    if (!is.null(note)) {
      msg <- sprintf("%s (%s)", msg, note);
    }
    msg <- sprintf("%s.", msg);
    throw(msg);
  }

  verbose && cat(verbose, "Number of tabular binary files located: ", 
                                                        length(pathnames));
  verbose && print(verbose, pathnames);


  verbose && enter(verbose, "Scanning for a valid file");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look for first possible valid match
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_along(pathnames)) {
    pathname <- pathnames[kk];
    verbose && enter(verbose, "File #", kk, " (", pathname, ")");

    # Create object
    res <- newInstance(static, pathname);

    # Correct number of units?
    if (!is.null(nbrOfUnits)) {
      if (nbrOfUnits(res) != nbrOfUnits) {
        res <- NULL;
      }
    }

    if (!is.null(res)) {
      verbose && cat(verbose, "Found a valid tabular binary file");
      verbose && exit(verbose);
      break;
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(res)) {
    queryStr <- paste(c(chipType, tags), collapse=",");
    throw("Failed to located a (valid) tabular binary file: ", queryStr);
  }

  verbose && print(verbose, res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Final validation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(nbrOfUnits)) {
    if (nbrOfUnits(res) != nbrOfUnits) {
      throw("The number of units in the loaded ", class(static)[1], " does not match the expected number: ", nbrOfUnits(res), " != ", nbrOfUnits);
    }
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: UnitNamesFile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("indexOfUnits", "AromaUnitTabularBinaryFile", function(this, names, ...) {
  # Map the unit names to the ones in the unit names file
  unf <- getUnitNamesFile(this);
  unitNames <- getUnitNames(unf);
  idxs <- match(names, unitNames);
  idxs;
}, protected=TRUE)



setMethodS3("allocateFromUnitNamesFile", "AromaUnitTabularBinaryFile", function(static, unf, ...) {
  # Argument 'unf':
  unf <- Arguments$getInstanceOf(unf, "UnitNamesFile");
  allocateFromUnitAnnotationDataFile(static, udf=unf, ...);
}, static=TRUE, protected=TRUE)


setMethodS3("allocateFromUnitAnnotationDataFile", "AromaUnitTabularBinaryFile", function(static, udf, path=NULL, tags=NULL, footer=list(), ...) {
  # Argument 'udf':
  udf <- Arguments$getInstanceOf(udf, "UnitAnnotationDataFile");

  # Output path
  if (is.null(path)) {
    chipTypeS <- getChipType(udf, fullname=FALSE);
    path <- file.path("annotationData", "chipTypes", chipTypeS);
  }
  path <- Arguments$getWritablePath(path);

  # Get platform
  platform <- getPlatform(udf);

  # Number of units
  nbrOfUnits <- nbrOfUnits(udf);

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(udf);

  # Exclude 'monocell' tags (AD HOC)
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);

  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create microarray tabular binary file
  allocate(static, filename=filename, path=path, nbrOfRows=nbrOfUnits, 
                                platform=platform, chipType=chipType, ...);
}, static=TRUE, protected=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: UnitNamesFile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2012-11-13
# o Explicitly declared "cached:.unf".
# 2011-11-19
# o Now byChipType() for AromaUnitTabularBinaryFile gives an error
#   message with more information on which file it failed to locate,
#   e.g. by specifying filename extension it looked for.
# 2011-03-28
# o BUG FIX: allocateFromUnitAnnotationDataFile() for 
#   AromaUnitTabularBinaryFile would include chip type tags in the
#   path, e.g. annotationData/chipTypes/GenomeWidesSNP_6,Full.
# 2011-02-28
# o STANDARDIZATION: Now the default output path for all 
#   allocateFromUnitAnnotationDataFile() is
#   annotationData/chipTypes/<chipType>/.  Before it was the same
#   directory as the original annotation data file, which may for
#   instance have been in a deeper subdirectory, or recently also
#   in a sibling root path.
# 2009-07-08
# o Added allocateFromUnitAnnotationDataFile() to AromaUnitTabularBinaryFile.
# 2009-05-12
# o Removed getUnitNamesFile() from AromaUnitTabularBinaryFile.
# 2009-02-10
# o Added byChipType() to AromaUnitTabularBinaryFile with option to 
#   validate/select by the number of units.
# 2008-07-09
# o Now AromaUnitTabularBinaryFile extends AromaMicroarrayTabularBinaryFile,
#   which contains a lot of the methods previously in this class.
# 2008-05-19
# o Added getPlatform().
# o Added platform-independent allocateFromUnitNamesFile() which now also
#   writes footer attribute 'platform'.
# 2008-02-13
# o Added and updated Rdoc comments.
# 2008-01-19
# o Now AromaUnitTabularBinaryFile gets the chip type from the file footer.
# o ROBUSTNESS: Now fromChipType() of AromaUnitTabularBinaryFile validates
#   that the number of units in the located file match the number of units
#   in the CDF located using the same search parameters.
# 2007-12-10
# o Currently a AromaUnitTabularBinaryFile (e.g. AromaUgpFile) does not
#   contain information about the "fullname" chip type, but only the basic
#   chip-type name, e.g. we cannot infer the full chip-type name from 
#   'GenomeWideSNP_5,Full,r2.ugp', but only 'GenomeWideSNP_5'. The fullname
#   should be the same as the full chip-type name of the CDF used to define
#   the the unit map, e.g. 'GenomeWideSNP_5,Full.CDF'.
#   We should add a header (or footer) field in the file format that 
#   indicates the full chip type.  
#   However, until that is done, the best we can do is to turn to the ad
#   hoc solution of scanning for the CDF with the longest matching fullname,
#   if both 'GenomeWideSNP_5,Full.CDF' and 'GenomeWideSNP_5.CDF' exists,
#   the we match the former to 'GenomeWideSNP_5,Full,r2.ugp'.  The fullname
#   chip type of the UGP is then full chip-type name of the CDF.  NOTE,
#   there is major drawback with this.  If the user deletes the "full" CDF,
#   the above approach would all of a sudden return a different full name!
# o Added clearCache().
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
