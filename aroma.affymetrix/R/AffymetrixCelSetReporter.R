###########################################################################/**
# @RdocClass AffymetrixCelSetReporter
#
# @title "The AffymetrixCelSetReporter class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFileSetReporter".}
#   \item{.setClass}{The name of the class of the input set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
#*/###########################################################################
setConstructorS3("AffymetrixCelSetReporter", function(..., .setClass="AffymetrixCelSet") {
  extend(AffymetrixFileSetReporter(..., .setClass=.setClass), "AffymetrixCelSetReporter"
  )
})



setMethodS3("as.character", "AffymetrixCelSetReporter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip type:", getChipType(this)));
  s <- c(s, paste("Number of arrays:", length(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the data set"
#
# \description{
#  @get "title" for which the reporter is generating output.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "AffymetrixCelSet".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "AffymetrixCelSetReporter", function(this, ...) {
  getFileSet(this);
})

setMethodS3("nbrOfArrays", "AffymetrixCelSetReporter", function(this, ...) {
  length(getDataSet(this));
}, protected=TRUE)


setMethodS3("getChipType", "AffymetrixCelSetReporter", function(this, ...) {
  cs <- getDataSet(this);
  cdf <- getCdf(cs);
  getChipType(cdf, ...);
}, protected=TRUE)


setMethodS3("getPath", "AffymetrixCelSetReporter", function(this, ...) {
  mainPath <- getMainPath(this);

  # Chip type
  chipType <- getChipType(this);

  # Get report set
  set <- getReportSet(this);

  # The full path
  path <- filePath(mainPath, chipType, set);
  path <- Arguments$getWritablePath(path);

  path;
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2007-03-19
# o Created from ArrayExplorer.R.
##############################################################################
