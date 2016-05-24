###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod extractExpressionSet
# @alias extractExpressionSet
#
# @title "Extracts an in-memory ExpressionSet object"
#
# \description{
#  @get "title" from a @see "ChipEffectSet" object.
#  Note that any modifications done to the extracted object will \emph{not}
#  be reflected in the original chip-effect set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Optional argument passed to @seemethod "extractMatrix".}
#   \item{logBase}{An @integer specifying the base to use when
#    log-transforming the signals.  If @NULL, the signals are not
#    transformed, but kept as is.}
#   \item{annotationPkg}{(optional) A @character string specifies a
#    Bioconductor (ChipDb, CDF environment or PDInfo) annotation package,
#    which then sets the 'annotation' slot of the returned object.
#    If \code{"ChipDb"}, \code{"cdf"} or \code{"PDInfo"} th corresponding
#    ChipDB, CDF environment or PDInfo package is inferred from
#    the CDF's chip type.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "Biobase::ExpressionSet-class" object.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractExpressionSet", "ChipEffectSet", function(this, ..., logBase=2, orderUnitsBy=c("asis", "lexicographic"), annotationPkg=NULL, verbose=FALSE) {
  .require <- require # To please R CMD check
  requireNamespace("Biobase") || throw("Package not loaded: Biobase");
  ExpressionSet <- Biobase::ExpressionSet


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getInteger(logBase, range=c(1,Inf));
  }

  # Argument 'orderUnitsBy':
  orderUnitsBy <- match.arg(orderUnitsBy);

  # Argument 'annotationPkg':
  if (!is.null(annotationPkg)) {
    annotationPkg <- Arguments$getCharacter(annotationPkg);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extract an ExpressionSet");

  verbose && print(verbose, this);
  verbose && cat(verbose, "Number of arrays: ", length(this));


  # Annotation string
  annotation <- character();
  if (!is.null(annotationPkg)) {
    verbose && enter(verbose, "Infer annotation string from annotation package");

    if (is.element(annotationPkg, c("ChipDb", "cdf", "PDInfo"))) {
      verbose && enter(verbose, "Infer annotation package name from CDF");
      cdfM <- getCdf(this);
      chipType <- getChipType(cdfM);
      chipType <- gsub(",monocell", "", chipType, fixed=TRUE);
      verbose && cat(verbose, "Chip type: ", chipType);

      if (is.element(annotationPkg, c("ChipDb"))) {
        annotationPkg <- affy::cleancdfname(chipType, addcdf=FALSE);
        annotationPkg <- sprintf("%s.db", annotationPkg);
      } else if (is.element(annotationPkg, c("cdf"))) {
        annotationPkg <- affy::cleancdfname(chipType);
      } else if (is.element(annotationPkg, c("PDInfo"))) {
        annotationPkg <- .cleanPlatformName(chipType);
      }

      verbose && cat(verbose, "Inferred annotation package name: ", annotationPkg);

      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Loading annotation package");
    verbose && cat(verbose, "Annotation package: ", annotationPkg);
    .require(annotationPkg, character.only=TRUE) || throw("Bioconductor annotation package not available: ", annotationPkg);
    verbose && exit(verbose);

    ns <- asNamespace(annotationPkg);

    # Infer type of annotation package, i.e. CDF of PDInfo package.
    if (regexpr("^pd[.]", annotationPkg) != -1L) {
      db <- ns[[annotationPkg]];
      if (is.null(db)) {
        throw("Unknown type of PDInfo annotation package: ", annotationPkg);
      }

      # Sanity checks
      dim <- .geometry(db);
      cdfM <- getCdf(this);
      cdf <- getMainCdf(cdfM);
      dim0 <- getDimension(cdf);
      if (!isTRUE(all.equal(dim, dim0, check.attributes=FALSE))) {
        throw(sprintf("The chip dimension of the requested annotation package ('%s') does not match the CDF: (%s) != (%s)", annotationPkg, paste(dim, collapse=", "), paste(dim0, collapse=", ")));
      }

      annotation <- .annotation(db);
      # Not needed anymore
      db <- dim <- dim0 <- NULL; # Not needed anymore
    } else if (regexpr("cdf$", annotationPkg) != -1L) {
      db <- ns[[annotationPkg]];
      if (is.null(db)) {
        throw("Unknown type of CDF annotation package: ", annotationPkg);
      }
      unitNames <- ls(envir=db);
      annotation <- gsub("[.]cdf$", "", annotationPkg);
      # Not needed anymore
      db <- unitNames <- NULL; # Not needed anymore
    } else if (regexpr("[.]db$", annotationPkg) != -1L) {
      db <- ns[[annotationPkg]];
      if (is.null(db) || !inherits(db, "ChipDb")) {
        throw("Unknown type of DB annotation package: ", annotationPkg);
      }
      annotation <- gsub("[.]db$", "", annotationPkg);
      # Not needed anymore
      db <- NULL; # Not needed anymore
    } else {
      throw("Unknown type of annotation package: ", annotationPkg);
    }

    annotation <- Arguments$getCharacter(annotation);
    verbose && cat(verbose, "Annotation string: ", annotation);

    # Not needed anymore
    ns <- NULL; # Not needed anymore
    verbose && exit(verbose);
  } # if (!is.null(annotationPkg))


  verbose && enter(verbose, "Reading data");
  Y <- extractMatrix(this, field="theta", drop=FALSE, returnUgcMap=TRUE, ..., verbose=less(verbose, 5));
  ugcMap <- attr(Y, "unitGroupCellMap");
  attr(Y, "unitGroupCellMap") <- NULL;
  verbose && str(verbose, Y);

  # "sdTheta" is actually standard errors, cf. rmaModelAffyPlm().
  # /HB 2014-04-27
  Yse <- extractMatrix(this, field="sdTheta", cells=ugcMap$cell, drop=FALSE, ..., returnUgcMap=FALSE, verbose=less(verbose, 5));
  verbose && str(verbose, Yse);

  verbose && str(verbose, ugcMap);
  verbose && exit(verbose);

  if (!is.null(logBase)) {
    verbose && enter(verbose, 'Log-transforming theta signals ("expressions")');
    verbose && cat(verbose, "Log base: ", logBase);
    Y <- log(Y, base=logBase);
    verbose && str(verbose, Y);
    # Note that affy::justRMA() returns se.exprs on the non-log scale
    # so don't log those here. /HB 2014-04-27
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Creating feature names");

  verbose && enter(verbose, "Reading (unit,group) names from CDF");
  cdf <- getCdf(this);
  verbose && print(verbose, cdf);

  ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap,
                                         verbose=less(verbose, 10));
  verbose && str(verbose, ugNames);
  verbose && exit(verbose);

  # Turn (<unit>,<group>) names into <unit>,<group> names
  names <- paste(ugNames$unitName, ugNames$groupName, sep=",");
  names <- gsub(",$", "", names);
  verbose && str(verbose, names);
  verbose && exit(verbose);

  rownames(Y) <- names;
  rownames(Yse) <- names;

  if (orderUnitsBy == "lexicographic") {
    o <- order(names);
    Y <- Y[o,,drop=FALSE];
    Yse <- Yse[o,,drop=FALSE];
    o <- NULL;
  }

  # Not needed anymore
  # Not needed anymore
  names <- ugNames <- ugcMap <- cdf <- NULL;

  verbose && enter(verbose, "Allocating ExpressionSet object");
  assayData <- new.env();
  assayData$exprs <- Y;
  assayData$se.exprs <- Yse;
  eset <- ExpressionSet(assayData=assayData, annotation=annotation);
  verbose && print(verbose, eset);

  # Not needed anymore
  Y <- NULL; # Not needed anymore
  verbose && exit(verbose);

  verbose && exit(verbose);

  eset;
}) # extractExpressionSet()


###########################################################################
# HISTORY:
# 2014-04-27 [HB]
# o Added argument 'orderUnitsBy' to extractExpressionSet() for
#   ChipEffectSet.
# o Now extractExpressionSet() for ChipEffectSet also returns std. errors.
# o Updated extractExpressionSet() to use the ExpressionSet constructor.
# 2013-04-27 [HB]
# o Added argument 'annotationPkg' to extractExpressionSet() for
#   ChipEffectSet, which (indirectly) sets the 'annotation' slot
#   of the returned ExpressionSet.
# 2011-07-14 [HB]
# o Added argument 'logBase' to extractExpressionSet().
# 2011-07-09 [HB]
# o Added extractExpressionSet() for ChipEffectSet.
# o Created.
###########################################################################
