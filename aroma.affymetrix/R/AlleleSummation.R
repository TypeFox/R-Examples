###########################################################################/**
# @RdocClass AlleleSummation
#
# @title "The AlleleSummation class"
#
# \description{
#  @classhierarchy
#
#  This class takes allele-specific chip effect estimates of a
#  SnpChipEffectSet and returns a CnChipEffectSet holding the summed
#  allele estimates.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "SnpChipEffectSet".}
#   \item{ignoreNAs}{If @TRUE, missing values are excluded when summing
#      the signals from the two alleles.}
#   \item{...}{Arguments passed to @see "UnitModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AlleleSummation", function(dataSet=NULL, ignoreNAs=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "SnpChipEffectSet");
  }

  extend(UnitModel(dataSet=dataSet, ...), "AlleleSummation",
    ignoreNAs = ignoreNAs,
    "cached:.outputSet" = NULL
  )
})


setMethodS3("getAsteriskTags", "AlleleSummation", function(this, collapse=NULL, ...) {
  # Returns 'U' (but allow for future extensions)
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "SA";

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getRootPath", "AlleleSummation", function(this, ...) {
  "plmData";
}, protected=TRUE)


setMethodS3("findUnitsTodo", "AlleleSummation", function(this, ...) {
  outSet <- getChipEffectSet(this);
  findUnitsTodo(outSet, ...);
})


###########################################################################/**
# @RdocMethod getChipEffectSet
#
# @title "Gets the set of chip effects for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "ChipEffectSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffectSet", "AlleleSummation", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # The output set
  outputSet <- this$.outputSet;
  if (!is.null(outputSet))
    return(outputSet);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect set. The input data set is empty.");

  verbose && enter(verbose, "Getting chip-effect set from data set");
  cdfM <- getCdf(ds);

  # Gets the ChipEffects Class object
  clazz <- getChipEffectSetClass(this);
  outputSet <- clazz$fromDataSet(dataSet=ds, path=getPath(this), cdf=cdfM,
                                                    verbose=less(verbose));
  setMergeStrands(outputSet, getMergeStrands(ds));
  setCombineAlleles(outputSet, TRUE);
  verbose && exit(verbose);

  # Store in cache
  this$.outputSet <- outputSet;

  outputSet;
})


setMethodS3("getChipEffectSetClass", "AlleleSummation", function(static, ...) {
  CnChipEffectSet;
}, static=TRUE, private=TRUE)


setMethodS3("process", "AlleleSummation", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Summing allele-specific estimates");
  outputSet <- getChipEffectSet(this);
  units <- findUnitsTodo(this, verbose=less(verbose, 5));

  if (length(units) == 0) {
    verbose && cat(verbose, "Already done.");
    verbose && exit(verbose);
    return(outputSet);
  }

  cdf <- getCdf(this);
  inputSet <- getDataSet(this);
  verbose && print(verbose, inputSet);

  verbose && enter(verbose, "Units to be updated");
  verbose && str(verbose, units);
  unitNames <- getUnitNames(cdf, units=units);
  verbose && str(verbose, unitNames);
  verbose && exit(verbose);

## OLD:
## snps <- indexOf(cdf, "SNP");
  types <- getUnitTypes(cdf, units=units);
  snps <- which(types == 2);
  # Not needed anymore
  types <- NULL;

  # WORKAROUND: Some of the units reported as SNPs, may actually be
  # non-SNPs.  Keep only those with two groups
  nbrOfGroups <- nbrOfGroupsPerUnit(cdf, units=snps);
  ok <- (nbrOfGroups %in% c(2,4));
  snps <- snps[ok];
  # Not needed anymore
  ok <- nbrOfGroups <- NULL;

  otherUnits <- setdiff(units, snps);
  verbose && cat(verbose, "Non-SNP units:");
  verbose && str(verbose, otherUnits);

  snpUgcMap <- otherUgcMap <- NULL;

  ignoreNAs <- this$ignoreNAs;

  nbrOfArrays <- length(inputSet);
  for (aa in seq_len(nbrOfArrays)) {
    inputFile <- inputSet[[aa]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", aa, getName(inputFile), nbrOfArrays));
    outputFile <- outputSet[[aa]];

    if (length(otherUnits) > 0) {
      verbose && enter(verbose, "Copying signals for non-SNP units");
      if (is.null(otherUgcMap)) {
        verbose && enter(verbose, "Getting (unit, group, cell) map for non-SNPs");
        otherUgcMap <- getUnitGroupCellMap(inputFile, units=otherUnits, verbose=less(verbose,5));
        verbose && exit(verbose);
      }
      cells <- otherUgcMap[,"cell"];
      if (length(cells) > 0) {
        data <- .readCel(getPathname(inputFile), indices=cells,
                        readIntensities=TRUE, readStdvs=TRUE, readPixels=TRUE);
        data <- as.data.frame(data[c("intensities", "stdvs", "pixels")]);
        verbose && str(verbose, data);
        data <- cbind(cell=cells, data);
        updateDataFlat(outputFile, data=data);
        # Not needed anymore
        data <- NULL;
      } else {
        verbose && cat(verbose, "Nothing to do: All units are SNP units.");
      }
      # Not needed anymore
      cells <- NULL;
      verbose && exit(verbose);
    }

    if (length(snps) > 0) {
      verbose && enter(verbose, "Combining allele signals for SNP units");
      if (is.null(snpUgcMap)) {
        verbose && enter(verbose, "Getting (unit, group, cell) map for SNPs");
        snpUgcMap <- getUnitGroupCellMap(inputFile, units=snps, verbose=less(verbose, 5));
        verbose && exit(verbose);
      }
      cells <- snpUgcMap[,"cell"];
      data <- .readCel(getPathname(inputFile), indices=cells,
                      readIntensities=TRUE, readStdvs=TRUE, readPixels=FALSE);
      yAB <- data[["intensities"]];
      verbose && cat(verbose, "(yA,yB) signals:");
      verbose && str(verbose, yAB);
      sdAB <- data[["stdvs"]];
      verbose && cat(verbose, "Standard deviations (yA,yB) signals:");
      verbose && str(verbose, sdAB);
      # (A,B,A,B,A,B,...)
      yAB <- matrix(yAB, nrow=2);
      sdAB <- matrix(sdAB, nrow=2);

      # Sum the alleles
      y <- sd <- rep(NA_real_, ncol(yAB));
      okAB <- !is.na(yAB);
      # (1) No missing data
      ok <- okAB[1,] & okAB[2,];
      y[ok] <- yAB[1,ok] + yAB[2,ok];
      sd[ok] <- sqrt(sdAB[1,ok]^2 + sdAB[2,ok]^2);
      if (ignoreNAs) {
        # (2a) Missing data in allele A
        ok <- !okAB[1,] & okAB[2,];
        y[ok] <- yAB[2,ok];
        sd[ok] <- sdAB[2,ok];
        # (2b) Missing data in allele B
        ok <- okAB[1,] & !okAB[2,];
        y[ok] <- yAB[1,ok];
        sd[ok] <- sdAB[1,ok];
      }
      # Not needed anymore
      yAB <- sdAB <- NULL;

      verbose && cat(verbose, "y=yA+yB signals:");
      verbose && str(verbose, y);
      verbose && cat(verbose, "sd=sqrt(sdA^2+sdB^2) signals:");
      verbose && str(verbose, sd);
      # Store signals in the cell for the A alleles:
      cells <- matrix(cells, nrow=2);
      cells <- cells[1,];

      data <- cbind(cell=cells, intensities=y, stdvs=sd);
      # Not needed anymore
      cells <- y <- sd <- NULL;

      updateDataFlat(outputFile, data=data);
      # Not needed anymore
      data <- NULL;
      verbose && exit(verbose);
    } # if (length(snps) > 0)

    verbose && exit(verbose);
  } # for (aa ...)

  verbose && exit(verbose);

  outputSet;
})



############################################################################
# HISTORY:
# 2008-12-18
# o BUG FIX: process() of AlleleSummation used verbose but had no such
#   argument causing it to use the global one.
# 2008-12-13
# o More verbose output.
# 2008-12-10
# o Now process() returns immediately if already done.
# o Now process() only processes non-summed units.
# o Added findUnitsTodo() to AlleleSummation.
# o Now AlleleSummation is also estimating the pooled standard deviation.
# o BUG FIX: AlleleSummation would not work for chip types containing
#   exclusively SNP units.  It expected some non-SNP units.
# 2008-02-20
# o Created.
############################################################################
