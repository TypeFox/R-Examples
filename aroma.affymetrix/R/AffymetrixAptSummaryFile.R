setConstructorS3("AffymetrixAptSummaryFile", function(..., sep=",", .verify=TRUE) {
  this <- extend(AffymetrixTabularFile(..., .verify=FALSE), "AffymetrixAptSummaryFile");

#  if (.verify)
#    verify(this, ...);
  this;
})

setMethodS3("translateColumnNames", "AffymetrixAptSummaryFile", function(this, names, ...) {
  # Convert 'FOO_BAR.CEL' to 'FOO_BAR'
  names <- gsub("[.](cel|CEL)$", "", names);

  names;
}, protected=TRUE)


setMethodS3("getArrayNames", "AffymetrixAptSummaryFile", function(this, ...) {
  names <- getColumnNames(this);
  names <- names[-1];
  names;
})

setMethodS3("nbrOfArrays", "AffymetrixAptSummaryFile", function(this, ...) {
  length(getArrayNames(this, ...));
})

setMethodS3("getReadArguments", "AffymetrixAptSummaryFile", function(this, fileHeader=NULL, colClasses=c("*"="double", "probeset_id"="character"), ..., verbose=FALSE) {
  args <- NextMethod("getReadArguments", colClasses=colClasses, verbose=verbose);
  args$quote <- "";
  args$na.strings <- "";
  
  args;
}, protected=TRUE)


setMethodS3("getProbesetIds", "AffymetrixAptSummaryFile", function(this, force=FALSE, ...) {
  probesetIds <- this$.probesetIds;
  if (force || is.null(probesetIds)) {
    probesetIds <- readProbesetIds(this, ...);
    this$.probesetIds <- probesetIds;
  }
  probesetIds;
})

setMethodS3("readProbesetIds", "AffymetrixAptSummaryFile", function(this, ...) {
  data <- readDataFrame(this, colClasses=c("probeset_id"="character"), ...);
  data <- data[,1,drop=TRUE];
  data;
}, protected=TRUE)


setMethodS3("readArrays", "AffymetrixAptSummaryFile", function(this, patterns, ...) {
  # Argument 'patterns':
  patterns <- Arguments$getCharacters(patterns);
  patterns <- unique(patterns);
  
  colClasses <- rep("double", length(patterns));
  names(colClasses) <- patterns;
  
  data <- readDataFrame(this, colClasses=colClasses, ...);
  data <- as.matrix(data);
  attr(data, "quantificationScale") <- getQuantificationScale(this);
  
  data;
})



setMethodS3("readHeader", "AffymetrixAptSummaryFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  header <- NextMethod("readHeader", verbose=verbose);

  # Extract APT parameters from header comments
  comments <- header$comments;
  names(comments) <- sprintf("commentRow%03d", seq_along(comments));
  isParam <- (regexpr("^#%", comments) != -1);
  header$comments <- comments[!isParam];
  params <- comments[isParam];
  params <- gsub("^#%", "", params);
  params <- strsplit(params, split="=", fixed=TRUE);
  keys <- sapply(params, FUN=function(param) {
    param[1];
  });
  params <- lapply(params, FUN=function(param) {
    paste(param[-1], collapse="=");
  });
  names(params) <- keys;

  header <- c(list(params=params), header);
  
  header;
}, protected=TRUE)


setMethodS3("getParameter", "AffymetrixAptSummaryFile", function(this, name, ...) {
  header <- getHeader(this, ...);
  value <- header$params[[name]];
  if (!is.null(value))
    value <- trim(value);
  value;
}, protected=TRUE)


setMethodS3("getTimeStamp", "AffymetrixAptSummaryFile", function(this, format="%a %b %d %H:%M:%S %Y", ...) {
  value <- getParameter(this, "affymetrix-algorithm-param-apt-time-str");
  if (!is.null(format)) {
    value <- strptime(value, format=format);
  }
  value;
}, protected=TRUE)
 

setMethodS3("getChipType", "AffymetrixAptSummaryFile", function(this, ...) {
  getParameter(this, "affymetrix-algorithm-param-apt-opt-chip-type");
})


setMethodS3("getQuantificationScale", "AffymetrixAptSummaryFile", function(this, ...) {
  getParameter(this, "affymetrix-algorithm-param-quantification-scale");
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-01-13
# o Created.
############################################################################
