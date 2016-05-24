###########################################################################/**
# @set "class=PSCBS"
# @RdocMethod writeSegments
#
# @title "Writes the table of segments to file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name, tags}{Name and optional tags part of the filename}.
#   \item{path}{The directory where the file will be written.}
#   \item{addHeader}{If @TRUE, header comments are written.}
#   \item{createdBy}{A header comment of whom created the file.}
#   \item{splitters}{If @TRUE, each chromosome is separated by a row
#     of missing values.}
#   \item{overwrite, skip}{If an output file already exists, these
#     arguments specifies what should happen.}
#   \item{...}{Additional arguments pass to \code{getSegments()}.}
# }
#
# \value{
#   Returns the pathname of the the file written.
# }
# 
# @author "HB"
#
# \seealso{
#   Utilizes @seemethod "getSegments".
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################  
setMethodS3("writeSegments", "PSCBS", function(fit, name=getSampleName(fit), tags=NULL, ext="tsv", path=NULL, addHeader=TRUE, createdBy=NULL, sep="\t", nbrOfDecimals=4L, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...) {
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

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);


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

  # Write to temporary file
  pathnameT <- pushTemporaryFile(pathname);


  sampleName <- getSampleName(fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getSegments(fit, ..., splitters=splitters);

  # Round of floating points
  if (!is.null(nbrOfDecimals)) {
    cols <- tolower(colnames(data));
    isInt <- (regexpr("chromosome|start|end|nbrofloci", cols) != -1);
    cols <- which(isInt);
    for (cc in cols) {
      values <- data[[cc]];
      if (is.double(values)) {
        values <- round(values, digits=0);
        data[[cc]] <- values;
      }
    } # for (key ...)

    cols <- tolower(colnames(data));
    isInt <- (regexpr("chromosome|start|end|nbrofloci", cols) != -1);
    isLog <- (regexpr("call", cols) != -1);
    isDbl <- (!isInt & !isLog);
    cols <- which(isDbl);
    for (kk in cols) {
      values <- data[[kk]];
      if (is.double(values)) {
        values <- round(values, digits=nbrOfDecimals);
        data[[kk]] <- values;
      }
    } # for (key ...) 
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (addHeader) {
#    sigmaDelta <- estimateStandardDeviation(fit, method="diff");
   sigmaDelta <- NA;
#    sigmaResiduals <- estimateStandardDeviation(fit, method="res");

    createdOn <- format(Sys.time(), format="%Y-%m-%d %H:%M:%S %Z");
    hdr <- c(
      name=name,
      tags=tags,
      fullname=fullname,
      segmentationMethod=sprintf("segment() of %s", attr(fit, "pkgDetails")),
      nbrOfLoci=nbrOfLoci(fit),
      nbrOfSegments=nbrOfSegments(fit),
      joinSegments=fit$params$joinSegments,
#      signalType=getSignalType(fit),
      sigmaDelta=sprintf("%.4f", sigmaDelta),
#      sigmaResiduals=sprintf("%.4f", sigmaResiduals),
      createdBy=createdBy,
      createdOn=createdOn,
      nbrOfDecimals=nbrOfDecimals,
      nbrOfColumns=ncol(data),
      columnNames=paste(colnames(data), collapse=", "),
      columnClasses=paste(sapply(data, FUN=function(x) class(x)[1]), collapse=", ")
    );
    bfr <- paste("# ", names(hdr), ": ", hdr, sep="");

    cat(file=pathnameT, bfr, sep="\n");
  } # if (addHeader)

  write.table(file=pathnameT, data, append=TRUE, quote=FALSE, sep=sep, 
                                          row.names=FALSE, col.names=TRUE);

  pathname <- popTemporaryFile(pathnameT);

  pathname;  
}) # writeSegments()



############################################################################
# HISTORY:
# 2011-12-03
# o Added writeSegments() for PSCBS.
############################################################################
