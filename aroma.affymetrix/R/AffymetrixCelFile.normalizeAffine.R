setMethodS3("transformAffine", "AffymetrixCelFile", function(this, outPath=file.path("transAffine", getChipType(this)), offset=0, scale=1, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  outPath <- Arguments$getWritablePathname(outPath);
  if (identical(getPath(this), outPath)) {
    throw("Cannot not transform data. Argument 'outPath' refers to the same path as the path of the data file to be transformed: ", outPath);
  }

  # Argument 'offset':
  offset <- Arguments$getDouble(offset);

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  cdf <- getCdf(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this));
  filename <- gsub("[.]cel$", ".CEL", filename);  # Only output upper case!
  pathname <- Arguments$getWritablePathname(filename, path=outPath,
                                         mustNotExist=(!overwrite && !skip));
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

  # Already shifted?
  if (skip && isFile(pathname)) {
    verbose && cat(verbose, "Transformed data file already exists: ", pathname);
    # CDF inheritance
    res <- fromFile(this, pathname);
    setCdf(res, cdf);
    return(res);
  }

  # Get probe signals
  x <- getData(this, fields="intensities", ..., verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subsetToUpdate <- identifyCells(cdf, probes=subsetToUpdate, types=typesToUpdate, verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift intensities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Transforming probe intensities by (offset,scale)=(%.1f,%.2f) ", offset, scale));
  x[subsetToUpdate] <- offset + scale*x[subsetToUpdate];
  # Not needed anymore
  subsetToUpdate <- NULL;
  verbose && exit(verbose);

  # Write normalized data to file
  verbose && enter(verbose, "Writing transformed probe signals");

  # Write to a temporary file (allow rename of existing one if forced)
  isFile <- (!skip && isFile(pathname));
  pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathnameT, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  .updateCel(pathnameT, intensities=x);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);

  # Return transformed data file object
  # CDF inheritance
  res <- fromFile(this, pathname);
  setCdf(res, cdf);
  return(res);
}, private=TRUE) # transformAffine()



############################################################################
# HISTORY:
# 2006-10-06
# o make sure cdf association is inherited
# 2006-08-25
# o Move to class AffymetrixCelFile and output is now CEL files only.
# 2006-07-28
# o Added transformAffine().
############################################################################
