###########################################################################/**
# @RdocClass Model
#
# @title "The Model class"
#
# \description{
#  @classhierarchy
#
#  This class is abstract and represents a generic model that applies
#  to a data set.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{The data set to which this model should be fitted.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
#   \item{.onUnknownArgs}{A @character string specifying what should occur
#      if there are unknown arguments in \code{...}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("Model", function(dataSet=NULL, tags="*", ..., .onUnknownArgs=c("error", "warning", "ignore")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AffymetrixCelSet");
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    if (is.element(.onUnknownArgs, c("error", "warning"))) {
      argsStr <- paste(names(args), collapse=", ");
      msg <- sprintf("Unknown arguments: %s", argsStr);
      if (.onUnknownArgs == "error") {
        throw(msg);
      } else if (.onUnknownArgs == "warning") {
        warning(msg);
      }
    }
  }

  this <- extend(Object(), c("Model", uses("ParametersInterface")),
    .dataSet = dataSet,
    .tags = NULL
  );

  # Interpret and append tags
  setTags(this, tags);

  this;
}, abstract=TRUE)



setMethodS3("as.character", "Model", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  ds <- getDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  s <- c(s, sprintf("Input tags: %s", getTags(ds, collapse=",")));
  s <- c(s, sprintf("Output tags: %s", getTags(this, collapse=",")));
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title".
#  By default, this is the string \code{"model"} appended by the capitalized
#  name of the model class, e.g. \code{"modelUnit"}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getRootPath", "Model", function(this, ...) {
  # Default root path
  paste("model", class(this)[1], sep="");
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the name of the input data set, unless
#  an alias is set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seemethod "getAlias"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "Model", function(this, ...) {
  name <- getAlias(this);

  if (is.null(name)) {
    ds <- getDataSet(this);
    name <- getName(ds);
  }

  name;
})


###########################################################################/**
# @RdocMethod getAlias
#
# @title "Gets the name alias for the model"
#
# \description{
#  @get "title", if set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string, or @NULL.
# }
#
# \seealso{
#   @seemethod "getName"
#   @seemethod "setAlias"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlias", "Model", function(this, ...) {
  this$.alias;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the name alias for the model"
#
# \description{
#  @get "title", if set.
# }
#
# @synopsis
#
# \arguments{
#  \item{alias}{A @character string, or @NULL.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seemethod "getAlias"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setAlias", "Model", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    # An alias must be a valid filename
    alias <- Arguments$getFilename(alias);

    # Assert that no commas are used.
    if (regexpr("[,]", alias) != -1) {
      throw("Aliases (names) must not contain commas: ", alias);
    }
  }

  this$.alias <- alias;

  invisible(this);
}, protected=TRUE)


setMethodS3("getAsteriskTags", "Model", function(this, collapse=NULL, ...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  tags <- name;

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which consists of the tags of the input data set followed
#  by an additional set of tags added by the model.
# }
#
# @synopsis
#
# \arguments{
#  \item{collapse}{A @character string used to concatenate the tags.
#     If @NULL, the tags are not concatenated.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# \seealso{
#   @seemethod "setTags"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "Model", function(this, collapse=NULL, ...) {
  # "Pass down" tags from the input data set
  ds <- getDataSet(this);
  inputTags <- getTags(ds);

  # Get class specific tags
  tags <- this$.tags;

  # In case this$.tags is not already split
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Expand asterisk tags
  if (any(tags == "*")) {
    tags[tags == "*"] <- getAsteriskTags(this, collapse=",");
  }

  # Combine input tags and local tags
  tags <- c(inputTags, tags);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    if (length(tags) > 0)
      tags <- unlist(strsplit(tags, split=","));
  }

  # No tags?
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})



###########################################################################/**
# @RdocMethod setTags
#
# @title "Sets the tags to be appended"
#
# \description{
#  @get "title" to the tags of the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{tags}{A @character @vector of tags.
#    The tags may also be passed as comma-separated strings.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seemethod "getTags"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setTags", "Model", function(this, tags=NULL, ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getTags(tags, collapse=NULL);
  }

  this$.tags <- tags;
})


###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output set"
#
# \description{
#  @get "title", which consists of the name with appended
#  comma-separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "Model", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of this model"
#
# \description{
#  @get "title" where the parameter files are located.
#  The directory is created, if missing.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \details{
#  If the path does not exist, it is created.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "Model", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType);
  path <- Arguments$getWritablePath(path);

  path;
})


###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the input data set for this model"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "AffymetrixCelSet".
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "Model", function(this, ...) {
  this$.dataSet;
})


###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "Model", function(this, ...) {
  getCdf(getDataSet(this));
}, private=TRUE)



###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments specific to any subclass.}
# }
#
# \value{
#  Returns an @integer @vector specifying what units where fitted.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "Model", abstract=TRUE);


setMethodS3("getLabel", "Model", function(this, ...) {
  label <- this$.label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
}, private=TRUE)

setMethodS3("setLabel", "Model", function(this, label, ...) {
  oldLabel <- this$.label;
  this$.label <- label;
  invisible(oldLabel);
}, private=TRUE)


############################################################################
# HISTORY:
# 2013-06-02
# o Added argument '.onUnknownArgs' to Model().
# 2008-09-03
# o CLEANUP: Removed 'parSet' argument of Model().
# 2007-10-11
# o Now getTags() of Model inserts "asterisk" tags.
# 2007-04-12
# o Added default getAsteriskTag().
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Added getAlias() and setAlias().
# 2007-01-14
# o Added a test for unknown arguments to constructor.  This was added
#   after long troubleshooting to find a call to MbeiPlm(mergeStrands=TRUE,
#   combineAlleles=TRUE) instead of MbeiCnPlm(...).
# 2007-01-06
# o Extracted from UnitModel.R.  Intended to cover all types of models,
#   e.g. RmaPlm, GladModel, PlasqModel etc.
# 2007-01-01
# o Created from former UnitGroupsModel with history as follows:
# 2006-11-19
# o Started to modify methods of this class to work similar to the
#   QuantileNormalizer and AllelicCrosstalkCalibrator classes.
# 2006-09-14
# o Not cloning the data set anymore.  Each model is responsible for
#   tranforming the data structure their way.  The advantage with this
#   approach is that we can cache read data in the data set object.
# 2006-08-28
# o Added getLabel(), which defaults to getName(), and setLabel().
# 2006-08-24
# o Added some Rdoc comments.
# 2006-08-17
# o Created.
############################################################################
