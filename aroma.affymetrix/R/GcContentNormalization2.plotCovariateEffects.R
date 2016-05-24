# Details:
# How is the Regional GC Correction algorithm enhancement implemented?
# The Regional GC Correction enhancement in the copy number algorithm utilizes information on GC content from the NetAffx NA26.1 annotation. For each marker, the annotation file reports the percent GC content in a half-megabase window. The CN algorithm bins the markers based on this GC content and normalizes the log2 ratios within each bin. Then the algorithm calculates an adjustment factor and applies this factor to each log2 ratio. This generates the GC-corrected log2 ratio for each marker. Please see the Copy Number Algorithm GC Waviness Correction white paper for complete details on the GC Correction algorithm enhancement.
# Reference:
# http://www.affymetrix.com/support/help/faqs/genotyping_console/copy_number_analysis/faq_5.jsp

setMethodS3("plotCovariateEffects", "GcContentNormalization2", function(this, arrays=NULL, units=NULL, ref="zero", ..., pch=".", lwd=2, xlim=NULL, ylim="auto", xlab="GC fraction", ylab="auto", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataSet <- getInputDataSet(this);
  # Argument 'arrays':
  if (!is.null(arrays)) {
    arrays <- Arguments$getIndices(arrays, max=length(dataSet));
  }

  if (is.null(ref)) {
  } else if (is.character(ref)) {
  } else if (inherits(ref, "SnpChipEffectFile")) {
  } else {
    throw("Unknown value of argument 'ref': ", class(ref)[1]);
  }


  unf <- getUnitNamesFile(dataSet);
  # Argument 'units':
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(unf));
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(unf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Plotting signals as a function of covariate(s)");
  verbose & cat(verbose, "Arrays:");
  verbose & str(verbose, arrays);

  X <- getCovariates(this, units=units, verbose=less(verbose,5));
  verbose & cat(verbose, "Covariates:");
  verbose & str(verbose, X);
  X <- X[,1];

  keep <- is.finite(X);
  X <- X[keep];
  units <- units[keep];
  # Not needed anymore
  keep <- NULL;
  # Sanity check
  stopifnot(length(units) == length(X));

  yR <- NULL;
  if (is.character(ref)) {
    if (ref == "median") {
      verbose && enter(verbose, "Calculating reference (average) file");
      ref <- getAverageFile(dataSet, verbose=less(verbose,10));
      verbose && exit(verbose);
      if (is.null(ylim)) {
        ylim <- c(-1,1)*3;
      }
    } else if (ref == "zero") {
      if (is.null(ylim)) {
        ylim <- c(0,16);
      }
    }
  }

  if (inherits(ref, "SnpChipEffectFile")) {
    verbose && enter(verbose, "Extracting reference signals");
    yR <- extractTotalAndFreqB(ref, units=units, drop=TRUE, verbose=less(verbose,50))[,"total"];
    yR <- log2(yR);
    verbose && str(verbose, yR);
    verbose && exit(verbose);
  }

  if (!is.null(arrays)) {
    dataSet <- extract(dataSet, arrays, onDuplicates="error");
  }

  if (is.null(xlim)) {
    xlim <- range(X, na.rm=TRUE);
  }

  if (identical(ylim, "auto")) {
    if (is.null(yR)) {
      ylim <- c(0,16);
    } else {
      ylim <- c(-3,3);
    }
  }

  if (identical(ylab, "auto")) {
    if (is.null(yR)) {
      ylab <- expression(log[2](theta));
    } else {
      ylab <- expression(log[2](theta/theta[R]));
    }
  }

  verbose && enter(verbose, "Plotting");
  nbrOfArrays <- length(dataSet);
  subplots(nbrOfArrays);
  par(mar=c(4,3,1,0.5)+0.1, mgp=c(2.0, 0.8, 0));
  for (cc in seq_len(nbrOfArrays)) {
    ce <- dataSet[[cc]];
    name <- getFullName(ce);
    name <- gsub(",chipEffects", "", name);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", cc, name, nbrOfArrays));
    x <- X;

    verbose && enter(verbose, "Extracting total signals");
    y <- extractTotalAndFreqB(ce, units=units, drop=TRUE, verbose=less(verbose,50))[,"total"];
    y <- log2(y);
    if (!is.null(yR)) {
      y <- y - yR;
    }
    verbose && exit(verbose);

    verbose && cat(verbose, "(x,y):");
    verbose && str(verbose, x);
    verbose && str(verbose, y);

    ok <- (is.finite(x) & is.finite(y));
    x <- x[ok];
    y <- y[ok];
    # Not needed anymore
    ok <- NULL;
    verbose && enter(verbose, "smoothScatter(x,y, ...):");
    verbose && str(verbose, x);
    verbose && str(verbose, y);
    smoothScatter(x, y, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
    verbose && exit(verbose);
    stext(side=3, pos=1, name);
    if (!is.null(yR)) {
      abline(h=0, col="#999999", lwd=lwd, lty=3);
    }
    fit <- smooth.spline(x,y);
    lines(fit, col="red", lwd=lwd);
    # Not needed anymore
    x <- y <- NULL;
    verbose && exit(verbose);
  } # for (cc ...)
  verbose && exit(verbose);
  # Not needed anymore
  yR <- NULL;

  verbose && exit(verbose);
}) # plotCovariateEffects()


############################################################################
# HISTORY:
# 2009-09-04
# o Now smoothScatter() is loaded via aroma.core.
# 2009-03-22
# o Created.
############################################################################
