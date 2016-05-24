###########################################################################/**
# @RdocClass UnitAnnotationDataFile
#
# @title "The UnitAnnotationDataFile interface class"
#
# \description{
#  @classhierarchy
#
#  A UnitAnnotationDataFile provides methods for querying certain types
#  of chip type annotation data by units.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.oo::Interface".}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("UnitAnnotationDataFile", function(...) {
  extend(Interface(), "UnitAnnotationDataFile");
})

setMethodS3("getChipType", "UnitAnnotationDataFile", function(...) {
  NextMethod("getChipType");
})

setMethodS3("getPlatform", "UnitAnnotationDataFile", function(...) {
  NextMethod("getPlatform");
})

setMethodS3("nbrOfUnits", "UnitAnnotationDataFile", function(...) {
  NextMethod("nbrOfUnits");
})

setMethodS3("getDefaultExtension", "UnitAnnotationDataFile", function(static, ...) {
  # Guess the filename extension from the class name, which might be wrong
  className <- class(static)[1];

  ext <- gsub("File$", "", className);
  ext <- strsplit(ext, split="", fixed=TRUE)[[1]];
  n <- length(ext);
  pos <- which(ext == toupper(ext));
  pos <- pos[length(pos)];
  ext <- ext[seq(from=pos, to=n)];
  ext <- paste(ext, collapse="");
  ext <- tolower(ext);

  attr(ext, "note") <- sprintf("this may not be the correct extension as it was guessed from the class name '%s'", className);

  ext;
}, static=TRUE, protected=TRUE);


setMethodS3("byChipType", "UnitAnnotationDataFile", function(static, chipType, tags=NULL, nbrOfUnits=NULL, ..., verbose=FALSE) {
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
                           firstOnly=FALSE, ..., verbose=less(verbose, 5));
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

  verbose && cat(verbose, "Number of ", class(static)[1], " located: ",
                                                        length(pathnames));
  verbose && print(verbose, pathnames);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Look for first possible valid match
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Scanning for a valid file");

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
      verbose && cat(verbose, "Found a valid ", class(static)[1]);
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



setMethodS3("getAromaUgpFile", "UnitAnnotationDataFile", function(this, ..., validate=FALSE, force=FALSE) {
  ugp <- this$.ugp;
  if (force || is.null(ugp)) {
    chipType <- getChipType(this, ...);
    ugp <- AromaUgpFile$byChipType(chipType, nbrOfUnits=nbrOfUnits(this), validate=validate);
    # Sanity check
    if (nbrOfUnits(ugp) != nbrOfUnits(this)) {
      throw("The number of units in located UGP file ('", getPathname(ugp), "') is not compatible with the data file ('", getPathname(this), "'): ", nbrOfUnits(ugp), " != ", nbrOfUnits(this));
    }
    this$.ugp <- ugp;
  }
  ugp;
})



setMethodS3("getAromaUflFile", "UnitAnnotationDataFile", function(this, ..., validate=FALSE, force=FALSE) {
  ufl <- this$.ufl;
  if (force || is.null(ufl)) {
    chipType <- getChipType(this, ...);
    ufl <- AromaUflFile$byChipType(chipType, nbrOfUnits=nbrOfUnits(this), validate=validate);
    # Sanity check
    if (nbrOfUnits(ufl) != nbrOfUnits(this)) {
      throw("The number of units in located UFL file ('", getPathname(ufl), "') is not compatible with the data file ('", getPathname(this), "'): ", nbrOfUnits(ufl), " != ", nbrOfUnits(this));
    }
    this$.ufl <- ufl;
  }
  ufl;
})



############################################################################
# HISTORY:
# 2011-11-19
# o Now byChipType() for UnitAnnotationDataFile gives an error
#   message with more information on which file it failed to locate,
#   e.g. by specifying filename extension it looked for.
# o Added default getDefaultExtension() for UnitAnnotationDataFile,
#   which guesses the filename extension from the class name.
# 2009-11-20
# o Now the "abstract" methods of the interface call NextMethod() instead,
#   which will give an error if nothing is defined.
# 2009-11-11
# o Added getAromaUflFile() to UnitAnnotationDataFile.
# 2009-07-08
# o Extracted methods from the UnitNamesFile interface class.
# o Created from UnitNamesFile.R.
############################################################################
