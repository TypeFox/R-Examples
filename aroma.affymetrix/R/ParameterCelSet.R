###########################################################################/**
# @RdocClass ParameterCelSet
#
# @title "The ParameterCelSet class"
#
# \description{
#  @classhierarchy
#
#  A ParameterCelSet object represents a set of @see "ParameterCelFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword "IO"
#*/###########################################################################
setConstructorS3("ParameterCelSet", function(...) {
  extend(AffymetrixCelSet(...), c("ParameterCelSet", uses("ParametersInterface")))
})


###########################################################################/**
# @RdocMethod extractMatrix
#
# @title "Extract data as a matrix for a set of arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{(The subset of units to be matched.
#     If @NULL, all units are considered.}
#   \item{...}{Passed to @see "base::subset" operating on the UGC map.}
#   \item{field}{The field to be extracted.}
#   \item{returnUgcMap}{If @TRUE, the (unit, group, cell) map is returned
#     as an attribute.}
#   \item{drop}{If @TRUE, singleton dimensions are dropped.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an JxK @double @matrix where J is the number of units,
#  and K is the number of arrays.
#  The names of the columns are the names of the arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractDataFrame".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMatrix", "ParameterCelSet", function(this, units=NULL, ..., field=c("intensities", "stdvs", "pixels"), returnUgcMap=FALSE, drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  ugcMap <- NULL;

  if (is.null(units)) {
    nunits <- nbrOfUnits(cdf);
  } else if (inherits(units, "UnitGroupCellMap")) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nunits <- length(units);
  }

  # Argument 'field':
  if (length(field) > 1)
    field <- field[1];

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Settings
  gcArrayFrequency <- getOption(aromaSettings, "memory/gcArrayFrequency", 10);

  verbose && enter(verbose, "Getting data for the array set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (unit, group, cell) map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }
  ugcMap <- subset(ugcMap, ...);

  if (nrow(ugcMap) == 0)
    throw("Nothing to return.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrayNames <- getNames(this);
  nbrOfArrays <- length(arrayNames);
  if (field %in% c("pixels")) {
    naValue <- as.integer(NA);
  } else {
    naValue <- as.double(NA);
  }
  df <- matrix(naValue, nrow=nrow(ugcMap), ncol=nbrOfArrays);
  colnames(df) <- arrayNames;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  for (aa in seq_len(nbrOfArrays)) {
    verbose && printf(verbose, "Array %d,\n", aa);
    cf <- this[[aa]];
    df[,aa] <- getDataFlat(cf, units=ugcMap, fields=field,
                                            verbose=less(verbose))[,field];
    if (aa %% gcArrayFrequency == 0) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }
  } # for (aa in ...)
  verbose && exit(verbose);

  # Drop singleton dimensions?
  if (drop) {
    df <- drop(df);
  }

  if (returnUgcMap)
    attr(df, "unitGroupCellMap") <- ugcMap;

  verbose && exit(verbose);

  df;
}) # extractMatrix()



###########################################################################/**
# @RdocMethod extractDataFrame
#
# @title "Extract data as a data.frame for a set of arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @seemethod "extractMatrix".}
#   \item{addNames}{If @TRUE, the first two columns contain the
#     unit names and the group names according the the CDF, otherwise
#     those two columns are not included.}
#   \item{addUgcMap}{If @TRUE, the columns following the unit and
#     group names contains the (unit, group, cell) index map.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a Jx(2+3+K) @data.frame where J is the number of units,
#  and K is the number of arrays.
#  The first two columns, if \code{addNames=TRUE}, contains the
#  unit names and the group names.
#  The next three columns contains the (unit, group, cell) index map.
#  The last K columns named by the arrays contain the data for the K arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractMatrix".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractDataFrame", "ParameterCelSet", function(this, addNames=FALSE, addUgcMap=TRUE, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array set");
  data <- extractMatrix(this, ..., returnUgcMap=TRUE,
                                                   verbose=less(verbose, 1));

  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  if (addUgcMap) {
    verbose && enter(verbose, "Merging UGC map and extracted data");
    ugcMap <- as.data.frame(ugcMap);
    data <- cbind(ugcMap, data);

    if (addNames) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  }

  if (addNames) {
    verbose && enter(verbose, "Appending unit and group names from CDF");
    cdf <- getCdf(this);
    verbose && cat(verbose, "CDF chip type: ",
                                        getChipType(cdf, fullname=TRUE));
    ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap,
                                              verbose=less(verbose, 10));
    # Not needed anymore
    cdf <- ugcMap <- NULL;
    verbose && cat(verbose, "(unit, group) names: ");
    verbose && str(verbose, ugNames);

    ugNames <- as.data.frame(ugNames);
    data <- cbind(ugNames, data);
    # Not needed anymore
    ugNames <- NULL;

    verbose && exit(verbose);
  }

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && exit(verbose);

  data;
}) # extractDataFrame()



############################################################################
# HISTORY:
# 2012-11-20
# o Added getParametersAsString() to ParameterCelSet.  Used to be in
#   direct subclasses.
# o Added getParameters() to ParameterCelSet, which calls ditto of
#   the first file.  Used to be in direct subclasses.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: extractMatrix().
# 2008-07-16
# o Added argument drop=FALSE to extractDataFrame().
# 2008-07-09
# o Added argument drop=FALSE to extractMatrix().
# 2008-02-28
# o Now argument 'units' also can be a UnitGroupCellMap.
# 2008-02-22
# o Created for the main purpose of putting extractMatrix() here, which
#   is currently common for ChipEffectSet and FirmaSet.
# o Moved extractDataFrame() from ChipEffectSet to ParameterCelSet.
# 2008-02-12
# o Added arguments 'addNames' and 'addUgcMap' to extractDataFrame().
# 2008-02-05
# o Created extractDataFrame() for ChipEffectSet.
# 2007-05-26
# o Updated the Rdocs for extractMatrix().
# 2007-03-04
# o Created extractMatrix().
############################################################################
