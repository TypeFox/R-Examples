setMethodS3("migrateTool", "CopyNumberSegmentationModel", function(static, what=c("addMissingAsteriskTag"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # In case the function was called using an instance.
  if (!inherits(static, "Class")) {
    static <- Class$forName(class(static)[1]);
  }

  if (what == "addMissingAsteriskTag") {
    throw("Hohoho! Don't run this one. It is only here for future use.");

    # Find all non-abstract subclasses
    classes <- c(getKnownSubclasses(static), getName(static));
    for (cc in seq_along(classes)) {
      className <- classes[cc];
      clazz <- Class$forName(className);
      # Ignore abstract classes
      if (isAbstract(clazz))
        classes[cc] <- NA;
    }
    classes <- na.omit(classes);


    # Scan the output directories for those classes
    for (cc in seq_along(classes)) {
      className <- classes[cc];
      clazz <- Class$forName(className);
      obj <- newInstance(clazz);
      aTag <- getAsteriskTags(obj)[1];
      pattern <- sprintf(",%s(,|$)", aTag);
      oldRootPath <- sprintf("%sData", tolower(aTag));

      verbose && enter(verbose, sprintf("Step #%d of %d - Scanning %s/ for datasets to be renamed", cc, length(classes), oldRootPath));

      dirs <- list.files(oldRootPath, full.names=TRUE);
      # Nothing to do?
      if (length(dirs) == 0) {
        verbose && exit(verbose);
        next;
      }

      # Consider only directories
      dirs <- dirs[sapply(dirs, FUN=isDirectory)];
      # Nothing to do?
      if (length(dirs) == 0) {
        verbose && exit(verbose);
        next;
      }

      # Identify data sets without the asterisk tag
      missingATag <- (regexpr(pattern, dirs) == -1);
      dirs <- dirs[missingATag];
      # Nothing to do?
      if (length(dirs) == 0) {
        verbose && exit(verbose);
        next;
      }

      verbose && enter(verbose, "Adding asterisk tag (", aTag, ") to ", length(dirs), " data sets");
      newDirs <- paste(dirs, aTag, sep=",");
      for (jj in seq_along(dirs)) {
        file.rename(dirs[jj], newDirs[jj]);
        verbose && printf(verbose, "Renamed '%s' to '%s' in '%s'.\n", basename(dirs[jj]), basename(newDirs[jj]), dirname(dirs[jj]));
      }
      verbose && exit(verbose);

      verbose && exit(verbose);
    } # for (cc ...)
  } # if (what == "addMissingAsteriskTag")
}, static=TRUE) # migrateTool()


##############################################################################
# HISTORY:
# 2007-09-05
# o Added static migrateTool() for CopyNumberSegmentationModel.
##############################################################################
