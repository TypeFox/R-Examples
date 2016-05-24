setConstructorS3("SampleAnnotationSet", function(...) {
  extend(GenericDataFileSet(...), "SampleAnnotationSet")
})


setMethodS3("findSAFs", "SampleAnnotationSet", function(static, path, pattern="[.](saf|SAF)$", ...) {
  # Search all paths to the root path
  pathnames <- list();
  lastPath <- NA;
  depth <- 10;
  while(depth > 0 && !is.null(path) && !identical(path, lastPath)) {
    lastPath <- path;
    pathnames0 <- list.files(path=path, pattern=pattern, full.names=TRUE);
    pathnames0 <- sort(pathnames0);
    pathnames <- c(pathnames, list(pathnames0));
#    path <- getParent(path);
    path <- dirname(path);
    depth <- depth - 1;
  }

  # Return from top to bottom
  pathnames <- rev(pathnames);

  pathnames <- unlist(pathnames, use.names=FALSE);

  pathnames;
}, static=TRUE, private=TRUE)


setMethodS3("byPathnames", "SampleAnnotationSet", function(static, pathnames, ..., fileClass="SampleAnnotationFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  pathnames <- sapply(pathnames, FUN=function(pathname) {
    Arguments$getReadablePathnames(pathname, mustExist=TRUE);
  });

  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  dfStatic <- Arguments$getInstanceOf(dfStatic, "SampleAnnotationFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Defining ", length(pathnames), " files");
  files <- list();
  for (kk in seq_along(pathnames)) {
    if (as.logical(verbose)) cat(kk, ", ", sep="");
    df <- newInstance(dfStatic, pathnames[kk]);
    files[[kk]] <- df;
    if (kk == 1) {
      # Update the static class instance.  The reason for this is
      # that if the second file cannot be instanciated with the same
      # class as the first one, then the files are incompatible.
      # Note that 'df' might be of a subclass of 'dfStatic'.
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }
  if (as.logical(verbose)) cat("\n");

  # Create the file set object
  set <- newInstance(static, files, ...);

  verbose && exit(verbose);

  set;
}, static=TRUE, private=TRUE) # byPathnames()


setMethodS3("fromPath", "SampleAnnotationSet", function(static, path, pattern="[.](saf|SAF)$", ..., fileClass="SampleAnnotationFile", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'pattern':
  if (!is.null(pattern))
    pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  dfStatic <- Arguments$getInstanceOf(dfStatic, "SampleAnnotationFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  pathnames <- findSAFs(static, path=path, pattern=pattern, ...);

  verbose && enter(verbose, "Defining ", length(pathnames), " files");
  files <- list();
  for (kk in seq_along(pathnames)) {
    if (as.logical(verbose)) cat(kk, ", ", sep="");
    df <- newInstance(dfStatic, pathnames[kk]);
    files[[kk]] <- df;
    if (kk == 1) {
      # Update the static class instance.  The reason for this is
      # that if the second file cannot be instanciated with the same
      # class as the first one, then the files are incompatible.
      # Note that 'df' might be of a subclass of 'dfStatic'.
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }
  if (as.logical(verbose)) cat("\n");
  verbose && exit(verbose);

  # Create the file set object
  set <- newInstance(static, files, ...);

  set;
}, static=TRUE) # fromPath()


setMethodS3("loadAll", "SampleAnnotationSet", function(static, paths="annotationData(|,.*)/", ..., merge=TRUE, reversePaths=TRUE, dropDuplicates=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paths':
  paths <- sapply(paths, FUN=Arguments$getRegularExpression);

  # Argument 'merge':
  merge <- Arguments$getLogical(merge);

  # Argument 'reversePaths':
  reversePaths <- Arguments$getLogical(reversePaths);

  # Argument 'dropDuplicates':
  dropDuplicates <- Arguments$getLogical(dropDuplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Loading all ", class(static)[1], ":s");

  verbose && enter(verbose, "Identifying all directories containing SAF files");
  pathnames <- findAnnotationData(set="samples", pattern="[.](saf|SAF)$",
                                firstOnly=FALSE, verbose=less(verbose,5));

  # Nothing to do?
  if (length(pathnames) == 0L) {
    verbose && cat(verbose, "No SAF files located.");
    verbose && exit(verbose);
    verbose && exit(verbose);
    return(list());
  }

  # Pathnames are now ordered according to aroma search conventions
  verbose && cat(verbose, "All SAF files located:");
  verbose && print(verbose, pathnames);

  paths <- dirname(pathnames);
  paths <- unique(paths);
  verbose && cat(verbose, "All directories with SAF files:");
  verbose && print(verbose, paths);
  verbose && exit(verbose);

  verbose && enter(verbose, "Loading ", class(static)[1], ":s");
  dsList <- lapply(paths, FUN=function(path) {
    fromPath(static, path=path, ..., verbose=verbose);
  })
  verbose && print(verbose, dsList);
  verbose && exit(verbose);

  if (merge) {
    # However, if we want to apply SAF files, we need to make sure
    # annotationData/ has higher priority than annotationData,<tags>/,
    # so we need to reverse the root paths order before merging.
    if (reversePaths) {
      verbose && enter(verbose, "Reversing order of ", class(static)[1], ":s");
      dsList <- rev(dsList);
      verbose && print(verbose, dsList);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Merging ", class(static)[1], ":s");
    dsList <- Reduce(append, dsList);
    verbose && print(verbose, dsList);
    verbose && exit(verbose);

    if (dropDuplicates) {
      verbose && enter(verbose, "Dropping duplicated files in ", class(static)[1], ":s");
      # (a) A duplicate must have the same file size as another file...
      ds <- dsList;

      # AD HOC: Undoing the above
      if (reversePaths) ds <- extract(ds, rev(seq_along(ds)));

      fileSizes <- sapply(ds, getFileSize);
      tt <- table(fileSizes);
      tt <- tt[tt > 1];
      dupFileSizes <- as.numeric(names(tt));
      dups <- (is.element(fileSizes, dupFileSizes));
      if (any(dups)) {
        # (b) ...and the same checksum
        dsT <- extract(ds, dups);
        checksumsT <- sapply(dsT, getChecksum);
        dupsT <- duplicated(checksumsT);

        if (any(dupsT)) {
          # Identify which to drop
          dups <- which(dups);
          dups <- dups[dupsT];

          verbose && printf(verbose, "Dropping %d files that are identical to other files:\n", length(dups));
          verbose && print(verbose, paste(getPathnames(ds)[dups]));
          ds <- extract(ds, -dups);

          # AD HOC: Redoing the above
          if (reversePaths) ds <- extract(ds, rev(seq_along(ds)));

          dsList <- ds;

          verbose && print(verbose, dsList);
        }
      }
      verbose && exit(verbose);
    }
  } # if (merge)

  verbose && exit(verbose);

  dsList;
}, static=TRUE, protected=TRUE) # loadAll()


############################################################################
# HISTORY:
# 2014-06-24
# o BUG FIX: SampleAnnotationSet$loadAll() would give an error if
#   annotationData/samples/ didn't exist or did not contain any SAF
#   files.
# 2011-03-03
# o Added internal static byPathnames().
# o Added static loadAll() for SampleAnnotationSet.
# 2008-05-09
# o Now SampleAnnotationSet inherits from GenericDataFileSet and no longer
#   from AffymetrixFileSet.
# 2007-03-06
# o Created.
############################################################################
