setConstructorS3("SampleAnnotationFile", function(...) {
  this <- extend(GenericDataFile(...), c("SampleAnnotationFile",
                                          uses("FileCacheKeyInterface")),
    "cached:.db" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getExtensionPattern", "SampleAnnotationFile", function(static, ...) {
  "[.](saf|SAF)$";
}, static=TRUE, protected=TRUE)


setMethodS3("fromPath", "SampleAnnotationFile", function(static, path, pattern=getExtensionPattern(static), ...) {
#  pathnames <- findSAFs(static, path=path, pattern=pattern, ...);
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, ...);
  if (length(pathnames) == 0)
    return(NULL);
  pathname <- pathnames[1];
  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)


setMethodS3("readDataFrame", "SampleAnnotationFile", function(this, rows=NULL, force=FALSE, ...) {
  db <- this$.db;
  if (force || is.null(db)) {
    pathname <- getPathname(this);

    # Read all non-commented lines
    bfr <- readLines(pathname);
    excl <- grep("^[ ]*#", bfr);
    if (length(excl) > 0)
      bfr <- bfr[-excl];

    # Parse these as a DCF
    con <- textConnection(bfr);
    on.exit(close(con));
    db <- read.dcf(con);
    db <- gsub("[\n\r]", "", db);
    # Not needed anymore
    bfr <- NULL;

    this$.db <- db;
  }

  colnames(db) <- toCamelCase(colnames(db));

  if (!is.null(rows))
    db <- db[rows,,drop=FALSE];

  db;
}, protected=TRUE)


setMethodS3("getPatterns", "SampleAnnotationFile", function(this, ...) {
  db <- readDataFrame(this, ...);

  # Get sample name pattern
  patterns <- sprintf("^%s.*$", db[,"name"]);
  patterns <- gsub("\\^\\^", "^", patterns);
  patterns <- gsub("\\$\\.\\*\\$", "$", patterns);

  patterns;
}, protected=TRUE)

setMethodS3("matchPatterns", "SampleAnnotationFile", function(this, names, trim=FALSE, ...) {
  # Scan vector of names for matching patterns
  patterns <- getPatterns(this, ...);
  res <- lapply(patterns, FUN=function(pattern) {
    idxs <- grep(pattern, names);
    names(idxs) <- names[idxs];
    idxs;
  });
  names(res) <- patterns; # In case length(res) == 1 /HB 2007-03-06

  if (trim) {
    keep <- (sapply(res, FUN=length) > 0);
    res <- res[keep];
  }

  res;
}, protected=TRUE)


setMethodS3("applyTo", "SampleAnnotationFile", function(this, names, FUN, ..., verbose=FALSE) {
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  allPatterns <- getPatterns(this, ..., verbose=verbose);

  res <- matchPatterns(this, names, trim=TRUE);
  # Nothing do to?
  if (length(res) == 0)
    return(invisible());

  verbose && print(verbose, res);

  patterns <- names(res);
  verbose && print(verbose, patterns);
  verbose && print(verbose, allPatterns);
  rows <- match(patterns, allPatterns);
  # Nothing do to?
  if (length(rows) == 0)
    return(invisible());

  verbose && print(verbose, rows);

  db <- readDataFrame(this, rows=rows);
  cc <- setdiff(colnames(db), "name");
  db <- db[,cc,drop=FALSE];

  # Nothing do to?
  if (nrow(db) == 0 || ncol(db) == 0)
    return(invisible());

  for (kk in seq_along(res)) {
    record <- db[kk,,drop=TRUE];

    # Nothing to do?
    if (all(is.na(record)))
      next;

    args <- list(
      appliesTo = res[[kk]]
    );
    args <- c(args, as.list(record));
    args <- c(args, list(...));
    do.call(FUN, args=args);
  }
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-05-09
# o Now SampleAnnotationFile inherits from GenericDataFile and no longer
#   from AffymetrixFile.
# 2008-04-14
# o Renamed readData() to readDataFrame() for SampleAnnotationFile.
# 2007-04-12
# o BUG FIX: readData() of SampleAnnotationFile would open a text connection
#   without closing it.
# 2007-03-13
# o getPatterns() and matchPatterns() now matches full names.
# 2007-03-06
# o Total make over.
# 2007-01-26
# o Created.
############################################################################
