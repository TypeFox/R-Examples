###########################################################################/**
# @RdocClass UnitModel
#
# @title "The UnitModel class"
#
# \description{
#  @classhierarchy
#
#  This class is abstract and represents a generic unit model, i.e.
#  a model that is applied to each unit separately.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet" to which this model should
#      be fitted.}
#   \item{probeModel}{A @character string specifying how PM and MM values
#      should be modelled.  By default only PM signals are used.}
#   \item{shift}{An optional amount the signals should be shifted
#      (translated) before fitting the model.}
#   \item{...}{Arguments passed to the constructor of @see "Model".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("UnitModel", function(dataSet=NULL, probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm"), shift=0, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AffymetrixCelSet");
  }

  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));


  extend(Model(dataSet=dataSet, ...), "UnitModel",
    probeModel = probeModel,
    shift = shift
  );
}, abstract=TRUE)



setMethodS3("getParameters", "UnitModel", function(this, ...) {
  params <- NextMethod("getParameters");
  params$probeModel <- this$probeModel;
  params$shift <- this$shift;
  params;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "UnitModel", function(this, collapse=NULL, ...) {
  # Returns 'U' (but allow for future extensions)
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "U";

  # Add class-specific tags
  shift <- as.integer(round(this$shift));
  if (shift != 0) {
    tags <- c(tags, sprintf("%+d", shift));
  }
  probeModel <- this$probeModel;
  if (probeModel != "pm") {
    tags <- c(tags, probeModel);
  }

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Gets the cell indices unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probesets).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to \code{getCellIndices()}
#     of the @see "AffymetrixCdfFile" of the input data set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure consisting of CDF cell indices.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCellIndices", "UnitModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Identify the type of probes to read
  stratifyBy <- switch(this$probeModel, "pm"="pm", "mm"="mm",
                       "pm-mm"="pmmm", "min1(pm-mm)"="pmmm", "pm+mm"="pmmm");

  # Get the CDF cell indices
  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  verbose && enter(verbose, "Identifying CDF cell indices");
  verbose && cat(verbose, "Stratify by: ", stratifyBy);
  cells <- getCellIndices(cdf, ..., stratifyBy=stratifyBy,
                                                      verbose=less(verbose));
  verbose && exit(verbose);

  cells;
}, private=TRUE)




###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
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
#  Returns an @integer @vector of unit indices.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "UnitModel", abstract=TRUE);


setMethodS3("getFitUnitFunction", "UnitModel", abstract=TRUE, private=TRUE);


setMethodS3("getFitSingleCellUnitFunction", "UnitModel", function(this, ...) {
  NULL;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-10-09
# o Added getFitSingleCellUnitFunction().
# 2008-09-05
# o TYPO: The error message for argument 'dataSet' in UnitModel() returned
#   multiple strings, one per class.
# 2008-09-03
# o CLEANUP: Moved the handling of 'probeModel' and 'shift'  to UnitModel
#   from ProbeLevelModel
# 2007-01-06
# o Removed never-used setup() method.
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
