setMethodS3("setAttributeXY", "AromaMicroarrayDataFile", function(this, value, ...) {
  # Argument 'value':
  if (is.null(value)) {
    # Nothing todo?
    return();
  }

  pattern <- "^(X*)(Y*)$";
  if (regexpr(pattern, value) == -1) {
    throw("The value of argument 'value' is unrecognized: ", value);
  }

  # Parse and count
  n23 <- gsub(pattern, "\\1", value);
  n24 <- gsub(pattern, "\\2", value);
  n23 <- nchar(n23);
  n24 <- nchar(n24);
  setAttributes(this, n23=n23, n24=n24);
}, protected=TRUE)

setMethodS3("getAttributeXY", "AromaMicroarrayDataFile", function(this, ...) {
  n23 <- getAttribute(this, "n23", 0);
  n24 <- getAttribute(this, "n24", 0);
  xyTag <- paste(c(rep("X", n23), rep("Y", n24)), collapse="");
  xyTag;
}, protected=TRUE)

setMethodS3("hasAttributeXY", "AromaMicroarrayDataFile", function(this, values, ...) {
   xyTag <- getAttributeXY(this);
   (xyTag %in% values);
}, protected=TRUE)


setMethodS3("setAttributesByTags", "AromaMicroarrayDataFile", function(this, tags=getTags(this), ...) {
  # Split tags
  tags <- Arguments$getTags(tags, collapse=NULL);

  newAttrs <- NextMethod("setAttributesByTags", tags=tags);

  # Parse XY, XX, XXX etc tags
  values <- grep("^X*Y*$", tags, value=TRUE);
  if (length(values) > 0) {
    newAttrs <- c(newAttrs, setAttributeXY(this, values));
  }

  # Parse tri<chromosome> tags
  values <- grep("^tri([1-9]|[0-9][0-9]|X|Y)$", tags, value=TRUE);
  if (length(values) > 0) {
    values <- gsub("^tri", "", values);
    chromosomes <- Arguments$getChromosomes(values);
    keys <- sprintf("n%02d", chromosomes);
    newAttrs <- c(newAttrs, lapply(keys, FUN=function(key) {
      setAttribute(this, key, 3);
    }));
  }

  # Return nothing
  invisible(newAttrs);
}, protected=TRUE)




setMethodS3("getPloidy", "AromaMicroarrayDataFile", function(this, chromosome, defaultValue=NA, ...) {
  # Argument 'chromosome':
  chromosome <- Arguments$getChromosome(chromosome);

  key <- sprintf("n%02d", chromosome);
  value <- getAttribute(this, key, defaultValue=defaultValue);
  value;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-03-05
# o Added setAttributesByTags().
# o Added setAttributeXY(), getAttributeXY(), and hasAttributeXY().
# 2007-02-12
# o Now getData() is using do.call() because it is faster. Unused arguments
#   are still ignored.
# 2007-02-04
# o Now getData() is call readCel() using doCall() so that unused arguments
#   in '...' are ignored.
# 2007-02-03
# o BUG FIX: getTimestamp() assumed a fix location in the CEL v3 header,
#   but that did not work for dChip exported CEL files.  Now, a US date
#   pattern is assumed and searched for.
# 2007-01-12 /KS
# o Moved image270() and writeSpatial() to AffymetrixCelFile.PLOT.R.
# 2006-12-18 /KS
# o Add "takeLog" argument (logical) to image270.  If true, take the log2
#   before plotting.  Can be more informative than natural scale.
# 2006-12-14
# o Removed getSampleName() which gives the same as getName().
# 2006-12-11
# o Now the timestamp is also reported for singel CEL files.
# o BUG FIX: getHeaderV3() would throw an error if there was an empty V3
#   header fields.  This was the reason why getTimestamp() gave an error
#   on some 100K chips.
# 2006-12-01
# o Added getTimestamp().
# 2006-11-28
# o Arguments 'force' and 'cache' has to be in readUnits() to avoid being
#   passed from calls of subclasses.
# 2006-10-23
# o Update default value for argument 'fields' in getData().
# 2006-10-22
# o In order to speed up fromFile(), the CEL header is not read anymore.
# 2006-10-06
# o make sure cdf association is inherited
# 2006-08-28
# o Renamed getFields() to getData() because getFields() is "reserved"
#   for use in the Object class.
# 2006-08-27
# o Added nbrOfCells() because it is so common.
# o Added createFrom() which utilizes new functions copyFile() and
#   clearData(). It is also no longer static. This is more generic and
#   cleaner.  The new clearData() does also not require the CDF file
#   (in case that should be missing).
# 2006-08-25
# o Renamed getIntensities() to getFields() which returns a data frame.
# o Added image270() and writeSpatial().
# o Added methods "[" and "[[" mapped to readUnits().
# 2006-08-24
# o Added the option to specify an 'cdf' object, making it possible to
#   override the default CDF file according to the CEL header.  It is
#   important that all methods queries the AffymetrixCdfFile object
#   from getCdf() and not the one through the CEL header!
# o Added most Rdoc comments.
# 2006-07-21
# o Added readUnits().
# 2006-07-05
# o BUG FIX/WORKAROUND: Currently the affxparser code crash R if the file
#   is not a valid CEL file.  The best we can do now is to test that the
#   filename has suffix *.CEL.
# 2006-05-30
# o Added fromFile().
# 2006-03-30
# o Updated according to affxparser.
# 2006-03-23
# o Moved all SNP related methods into the new class AffymetrixSnpCelFile.
# 2006-03-18
# o Made probe indices one-based.
# 2006-03-04
# o Added support for remapping in readIntensities().  This is currently
#   not used for CEL files (only APD files), but was added for the future.
# 2006-03-02
# o Created.
############################################################################
