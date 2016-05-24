# @author "HB"
#
# \references{
#   [1] http://wiki.transvar.org/confluence/display/igbman/File+Formats
# }
setMethodS3("writeSgr", "AffymetrixCelSet", function(this, units=NULL, ..., tags=getTags(this), verbose=FALSE, fileExtension="sgr", fileSep="\t", nbrOfFigures=3) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  chrText <- paste(c(1:22, "X", "Y", "M"), sep="");
  names(chrText) <- 1:25;

  nbrOfArrays <- length(this);
  cdf <- getCdf(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading indices");
  indices <- getCellIndices(cdf, units=units, stratifyBy="pm");
  indices <- unlist(indices, use.names=FALSE);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read probe genomic location
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Gathering probe genomic locations");
  acp <- AromaCellPositionFile$byChipType(getChipType(cdf));

  ch <- paste("chr", chrText[ as.character(acp[indices,1,drop=TRUE]) ], sep="");
  pos <- acp[indices,2,drop=TRUE];
  verbose && exit(verbose);

  for (ii in seq_len(nbrOfArrays)) {
    cf <- this[[ii]];
    sampleName <- getName(cf);

    verbose && enter(verbose, sprintf("Gathering and writing data for ", sampleName));

    fullname <- paste(c(sampleName, tags), collapse=",");
    filename <- sprintf("%s.%s", fullname, fileExtension);

    data <- extractMatrix(cf, cells=indices, verbose=verbose);
    data <- log2(data);
    data <- round(data, digits=nbrOfFigures);
    data <- cbind(ch, pos, data);
    write.table(data, filename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=fileSep);

    verbose && exit(verbose);
  } # for (ii ...)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
})
