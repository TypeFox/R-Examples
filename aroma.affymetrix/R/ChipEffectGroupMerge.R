###########################################################################/**
# @RdocClass ChipEffectGroupMerge
#
# @title "The ChipEffectGroupMerge class"
#
# \description{
#  @classhierarchy
#
#  This class represents a method that merges chip effects across groups
#  unit by unit.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "ChipEffectSet".}
#   \item{fcn}{A @function used to "merge" groups.}
#   \item{...}{Additional arguments passed to the constructor of
#     @see "ChipEffectTransform".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword "internal"
#*/###########################################################################
setConstructorS3("ChipEffectGroupMerge", function(dataSet=NULL, fcn=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "ChipEffectSet");
  }

  this <- extend(ChipEffectTransform(dataSet, ...), "ChipEffectGroupMerge",
    .mergeFcn = NULL
  )

  setMergeFunction(this, fcn);

  this;
})


setMethodS3("getMergeFunction", "ChipEffectGroupMerge", function(this, ...) {
  this$.mergeFcn;
})

setMethodS3("setMergeFunction", "ChipEffectGroupMerge", function(this, fcn=NULL, ...) {
  # Argument 'fcn':
  if (is.null(fcn)) {
  } else if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", class(fcn)[1]);
  }

  if (!is.null(fcn)) {
    # Test the merge function
    for (size in 1:10) {
      y <- matrix(1000+1:(size*4), nrow=size);
      yOut <- fcn(y);
      if (!identical(dim(yOut), dim(y))) {
        throw("Function 'fcn' must not change the dimension of the data: ",
                                    paste(dim(yOut), collapse="x"), " != ",
                                              paste(dim(y), collapse="x"));
      }
    }
  }

  this$.mergeFcn <- fcn;
})


setMethodS3("getParameters", "ChipEffectGroupMerge", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    fcn = getMergeFunction(this)
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


setMethodS3("getCdf", "ChipEffectGroupMerge", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})



###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "ChipEffectGroupMerge", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Merging unit-group data across arrays");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Nothing to do?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mergeFcn <- getMergeFunction(this);
  if (is.null(mergeFcn)) {
    throw("Nothing to do. There is no merge function: NULL");
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already merged");
    verbose && exit(verbose);
    outputSet <- getOutputDataSet(this);
    return(invisible(outputSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ces <- getInputDataSet(this);

  cdf <- getCdf(this);

  # Get (and create) the output path
  path <- getPath(this);
  path <- Arguments$getWritablePath(path);

  # Fields to be merged
  fields <- c("theta", "sdTheta");

  overwrite <- force;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fl <- NULL;
  targetFcn <- NULL;
  map <- NULL;
  nbrOfArrays <- length(ces);
  res <- vector("list", nbrOfArrays);
  for (kk in seq_len(nbrOfArrays)) {
    ce <- ces[[kk]];
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')", kk, nbrOfArrays, getName(ce)));
#    verbose && print(verbose, ce);

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already merged. Skipping.");
      ceM <- fromFile(ce, pathname);

      # CDF inheritance
      setCdf(ceM, cdf);

      res[[kk]] <- ceM;
      verbose && exit(verbose);
      next;
    }

    # Get unit-to-cell (for optimized reading)?
    if (is.null(map)) {
      # Only loaded if really needed.
      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
      map <- getUnitGroupCellMap(ce, verbose=less(verbose));
      verbose && str(verbose, map);
      verbose && exit(verbose);

      # Get all possible unit sizes (number of groups per unit)
      uSizes <- sort(unique(map[,"group"]));
      verbose && cat(verbose, "Unit sizes (number of groups) identified:");
      verbose && print(verbose, uSizes);
    }

    # Extract signals for all units
    verbose && enter(verbose, "Getting signals");
    data <- getDataFlat(ce, units=map, fields=fields, verbose=less(verbose));
    verbose && str(verbose, data);
    verbose && exit(verbose);

    # For each unit size
    excl <- c();
    for (size in rev(uSizes)) {
      verbose && enter(verbose, "Unit size ", size);
      # Identify the units of that size
      idxs <- which(data[,"group"] == size);
      unitsS <- data[idxs, "unit"];
      # Exclude already used units
      unitsS <- setdiff(unitsS, excl);
      # For next round, do not use the same units
      excl <- c(excl, unitsS);

      # Get the subset of the data for such units
      idxs <- which(data[,"unit"] %in% unitsS);
      # Not needed anymore
      unitsS <- NULL;

      for (field in fields) {
        # Extract signals as a matrix where each column is one unit
        y <- data[idxs, field];
        y <- matrix(y, nrow=size);
        verbose && str(verbose, y);
        y <- mergeFcn(y);
        verbose && str(verbose, y);

        # Update data table
        data[idxs, field] <- as.vector(y);

        # Not needed anymore
        # Not needed anymore
        y <- NULL;
      }

      # Not needed anymore
      # Not needed anymore
      idxs <- NULL;

      verbose && exit(verbose);
    } # for (size in ...)
      # Not needed anymore
    # Not needed anymore
    excl <- NULL;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Copy CEL file and update the copy
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing merged data");
    verbose && cat(verbose, "Pathname: ", pathname);
    tryCatch({
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      ceOut <- createFrom(ce, filename=pathname, path=NULL,
                      methods="create", clear=TRUE, verbose=less(verbose));
      verbose && exit(verbose);

      verbose && print(verbose, ceOut);
      verbose && enter(verbose, "Writing merged data");
      updateDataFlat(ceOut, data=data, verbose=less(verbose));
      verbose && exit(verbose);
    }, error = function(cond) {
      verbose && enter(verbose, "Error detected while saving, removing data file again");
      verbose && print(verbose, cond);
      verbose && cat(verbose, "Pathname: ", pathname);
      verbose && exit(verbose);
      throw(cond$message);
    })
    verbose && exit(verbose);

    # Not needed anymore
    # Not needed anymore
    ce <- data <- NULL;

    # CDF inheritance
    setCdf(ceOut, cdf);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    res[[kk]] <- ceOut;

    verbose && exit(verbose);
  } # for (kk in ...)

  # Create the output set by re-reading it from files
  outputSet <- byPath(ces, path, cdf=cdf, ...);

  # Update the output data set
  this$outputSet <- outputSet;

  verbose && exit(verbose);

  outputSet;
})

############################################################################
# HISTORY:
# 2007-02-21
# o BUG FIX: When identifying units of certain sizes, previously used units
#   where not excluded resulting in groups =K, >=K-1, ..., >=1 rather than
#   =K, =K-1, ..., =1.
# 2007-02-20
# o Created.
############################################################################
