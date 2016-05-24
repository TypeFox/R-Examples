setMethodS3("extractWIG", "AbstractCBS", function(fit, signal, transform=NULL, nbrOfDecimals=4L, label=toupper(signal), graphType=c("bar", "points", "line"), viewLimits=NULL, colors=c(negative="231,41,138", positive="117,112,179"), ...) {
  # Argument 'graphType':
  graphType <- match.arg(graphType)

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);

  data <- getSegments(fit, splitter=FALSE)
  fields <- c("chromosome", "start", "end", "mean")
  if (!is.null(signal)) {
    fields[-1] <- sprintf("%s%s", signal, capitalize(fields[-1]))
  }
  data <- data[,fields]
  colnames(data) <- c("chromosome", "start", "end", "mean")
  data$chromosome <- sprintf("chr%d", data$chromosome)

  ## Round / truncate
  for (ff in c("start", "end")) {
    data[[ff]] <- as.integer(round(data[[ff]], digits=0L))
  }

  # Transform mean levels?
  if (!is.null(transform)) {
    data[["mean"]] <- transform(data[["mean"]])
  }
  
  # Round mean levels
  if (!is.null(nbrOfDecimals)) {
    data[["mean"]] <- round(data[["mean"]], digits=nbrOfDecimals);
  }

  # Drop segments with missing values
  data <- na.omit(data)

  ## Track information
  track <- list(
    type="wiggle_0",
    name=sampleName(fit),
    description=sprintf("Data type: %s", class(fit)),
    graphType=graphType,
    visibility="full",
    maxHeightPixels="128:96:64",
    yLineOnOff="on",
    autoScale="true"
  )
  if (is.na(track$name)) track$name <- "Unknown sample"
  if (!is.null(signal)) track$name <- sprintf("%s [%s]", track$name, label)

  if (!is.null(viewLimits)) {
    track$viewLimits <- sprintf("%g:%g", viewLimits[1], viewLimits[2])
  }

  if (!is.null(colors)) {
    if (!is.null(names(colors))) colors <- colors[c("negative", "positive")]
    track$color <- colors[["negative"]]
    track$altColor <- colors[["positive"]]
  }

  attr(data, "track") <- track

  data
}, protected=TRUE)



setMethodS3("extractWIG", "CBS", function(fit, ..., colors=c(negative="231,41,138", positive="117,112,179")) {
  NextMethod("extractWIG", signal=NULL, colors=colors)
}, protected=TRUE)


setMethodS3("extractWIG", "PSCBS", function(fit, signal=c("tcn", "dh"), ..., colors=c(negative="231,41,138", positive="117,112,179")) {
  signal <- match.arg(signal)
  NextMethod("extractWIG", signal=signal, colors=colors)
}, protected=TRUE)


# \references{
#  [1] Wiggle Track Format (WIG), UCSC Genome Browser
#      \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}
# }
setMethodS3("writeWIG", "AbstractCBS", function(fit, name=getSampleName(fit), tags=NULL, ext="wig", path=NULL, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name' and 'tags':
  name <- Arguments$getCharacter(name);
  tags <- Arguments$getCharacters(tags);

  # Argument 'ext':
  ext <- Arguments$getCharacter(ext);

  # Arguments 'path':
  path <- Arguments$getWritablePath(path);


  fullname <- paste(c(name, tags), collapse=",");
  filename <- sprintf("%s.%s", fullname, ext);
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # File already exists?
  if (isFile(pathname)) {
    # Skip?
    if (skip) {
      return(pathname);
    }

    # Overwrite!
    file.remove(pathname);
  }

  ## Write file (atomically)
  pathnameT <- pushTemporaryFile(pathname)

  bed <- extractWIG(fit, ...)

  ## Generate 'track' definition line
  track <- attr(bed, "track")
  track <- lapply(track, FUN=function(value) {
    if (is.character(value)) value <- dQuote(value)
    value
  })
  track <- unlist(track, use.names=TRUE)
  track <- sprintf("%s=%s", names(track), track)
  track <- paste(track, collapse=" ")
  track <- sprintf("track %s", track)


  cat(track, "\n", sep="", file=pathnameT)
  write.table(bed, file=pathnameT,
              col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE,
              append=TRUE)

  pathname <- popTemporaryFile(pathnameT)

  pathname
})


############################################################################
# HISTORY:
# 2015-09-08
# o Added extractWIG() and writeWIG() for CBS objects.
############################################################################
