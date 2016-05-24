###########################################################################/**
# @RdocClass AromaCellTabularBinaryFile
#
# @title "The AromaCellTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaCellTabularBinaryFile is an @see "AromaTabularBinaryFile" with
#  the constraint that the rows map one-to-one to the cells (features)
#  of a microarray.
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
# \seealso{
#   @see "AromaUnitTabularBinaryFile".
# }
#*/###########################################################################
setConstructorS3("AromaCellTabularBinaryFile", function(...) {
  extend(AromaMicroarrayTabularBinaryFile(...), "AromaCellTabularBinaryFile");
})


setMethodS3("nbrOfCells", "AromaCellTabularBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("byChipType", "AromaCellTabularBinaryFile",function(static, chipType, tags=NULL, nbrOfCells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'nbrOfCells':
  if (!is.null(nbrOfCells)) {
    nbrOfCells <- Arguments$getInteger(nbrOfCells, range=c(0,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating ", class(static)[1]);
  pathname <- findByChipType(static, chipType=chipType, tags=tags,
                                     firstOnly=TRUE, ...);
  if (is.null(pathname)) {
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

  verbose && cat(verbose, "Located file: ", pathname);
  res <- newInstance(static, pathname);
  verbose && print(verbose, res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validation?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(nbrOfCells)) {
    if (nbrOfCells(res) != nbrOfCells) {
      throw("The number of cells in the loaded ", class(static)[1], " does not match the expected number: ", nbrOfCells(res), " != ", nbrOfCells);
    }
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)



############################################################################
# HISTORY:
# 2014-06-28
# o Added byChipType() for AromaCellTabularBinaryFile.
# 2008-07-09
# o Created from AromaUnitTabularBinaryFile.R.
############################################################################
