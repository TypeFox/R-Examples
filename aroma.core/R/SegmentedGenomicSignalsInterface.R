###########################################################################/**
# @RdocClass SegmentedGenomicSignalsInterface
#
# @title "The SegmentedGenomicSignalsInterface class interface"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("SegmentedGenomicSignalsInterface", function(...) {
  extend(Interface(), "SegmentedGenomicSignalsInterface");
})

setMethodS3("setStates", "SegmentedGenomicSignalsInterface", function(this, states=NULL, ...) {
  # Argument 'states':
  if (!is.null(states)) {
    if (is.function(states)) {
    } else {
      throw("Argument 'states' must be a function: ", mode(states));
    }
  }

  if (inherits(this, "RichDataFrame")) {
    fcn <- function(data, ...) {
      x <- data$x;
      states(x, ...);
    }
    this[["state"]] <- fcn;
  } else {
    this <- setBasicField(this, ".states", states);
  }

  invisible(this);
})


setMethodS3("getStates", "SegmentedGenomicSignalsInterface", function(this, x=NULL, ...) {
  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getNumerics(x, disallow=NULL);
    nbrOfLoci <- length(x);
  }

  if (inherits(this, "RichDataFrame")) {
    if (is.null(x)) {
      states <- this$state;
      nbrOfLoci <- length(states);
    } else {
      tmp <- newInstance(this, nrow=length(x));
      tmp$x <- x;
      states <- tmp$state;
      # Not needed anymore
      tmp <- NULL;
    }
  } else {
    if (is.null(x)) {
      x <- getPositions(this);
    }
    nbrOfLoci <- length(x);

    states <- getBasicField(this, ".states");

    if (is.function(states)) {
      fcn <- states;
      chromosome <- getChromosome(this);
      name <- getName(this);
      states <- fcn(x, chromosome=chromosome, name=name, ...);
      storage.mode(states) <- "integer";
    }
  }

  # Sanity check
  stopifnot(length(states) == nbrOfLoci);

  states;
})


setMethodS3("getUniqueStates", "SegmentedGenomicSignalsInterface", function(this, na.rm=FALSE, ...) {
  states <- getStates(this, ...);
  states <- unique(states);
  na.last <- if (na.rm) NA else TRUE;
  states <- sort(states, na.last=na.last);
  states;
})


setMethodS3("findChangePointsByState", "SegmentedGenomicSignalsInterface", function(this, na.rm=FALSE, ends=FALSE, ...) {
  x <- this$x;
  s <- getStates(this, x=x);
  nas <- is.na(s);
  if (na.rm) {
    keep <- which(!nas);
    s <- s[keep];
    x <- x[keep];
  } else {
    s[nas] <- +Inf;
  }
  ds <- diff(s);
  idxs <- which(ds != 0);
  xCP <- (x[idxs] + x[idxs+1L]) / 2;
  if (ends) {
    xCP <- c(min(x, na.rm=TRUE), xCP, max(x, na.rm=TRUE));
  }
  xCP;
}) # findChangePointsByState()


setMethodS3("as.data.frame", "SegmentedGenomicSignalsInterface", function(x, ..., virtual=TRUE) {
  # To please R CMD check
  this <- x;

  if (inherits(this, "RichDataFrame")) {
    df <- NextMethod("as.data.frame", virtual=virtual);
  } else {
    df <- NextMethod("as.data.frame");
    if (virtual) {
      df$state <- getStates(this, x=df$x);
    }
  }

  df;
})

setMethodS3("getVirtualField", "SegmentedGenomicSignalsInterface", function(this, key, ...) {
  # Argument 'key':
  key <- Arguments$getCharacter(key);

  if (inherits(this, "RichDataFrame")) {
    value <- this[[key]];
  } else {
    fields <- getVirtualLocusFields(this, ...);
    stopifnot(is.element(key, fields));

    if (key == "state") {
      value <- getStates(this, ...);
    } else {
      value <- NextMethod("getVirtualField", key=key);
    }
  }

  value;
}, protected=TRUE)

setMethodS3("getVirtualLocusFields", "SegmentedGenomicSignalsInterface", function(this, ...) {
  if (inherits(this, "RichDataFrame")) {
    fields <- getVirtualColumnNames(this, ...);
  } else {
    fields <- NextMethod("getVirtualLocusFields");
    fields <- c(fields, "state");
    fields <- unique(fields);
  }
  fields;
}, protected=TRUE)


setMethodS3("extractSubsetByState", "SegmentedGenomicSignalsInterface", function(this, states, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'states':
  states <- Arguments$getVector(states);


  # Identify loci that have the requested states
  signalStates <- getStates(this);

  # Sanity check
  uniqueStates <- unique(signalStates);
  if (length(uniqueStates) > 0.1*length(signalStates)) {
    throw(sprintf("Detected too many states: %s [%d]", hpaste(signalStates), length(uniqueStates)));
  }

  # Subset
  keep <- is.element(signalStates, states);
  keep <- which(keep);

  # Extract this subset
  extractSubset(this, keep, ...);
})



setMethodS3("kernelSmoothingByState", "SegmentedGenomicSignalsInterface", function(this, xOut=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'xOut':
  if (!is.null(xOut)) {
    xOut <- Arguments$getNumerics(xOut);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  x <- getPositions(this);

  if (is.null(xOut)) {
    xOut <- x;
  }
  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  naValue <- as.double(NA);
  yOut <- rep(naValue, length(xOut));

  y <- getSignals(this);
  states <- getStates(this);

  uniqueStates <- unique(states);
  uniqueStates <- sort(uniqueStates, na.last=TRUE);
  verbose && cat(verbose, "Unique states:");
  verbose && str(verbose, uniqueStates);

  # Identify states of target loci
  statesOut <- getStates(this, x=xOut);

  for (ss in seq_along(uniqueStates)) {
    state <- uniqueStates[ss];
    verbose && enter(verbose, sprintf("State #%d ('%d') of %d",
                                      ss, state, length(uniqueStates)));

    # Identifying loci with this state
    if (is.na(state)) {
      keep <- is.na(states);
    } else {
      keep <- (states == state);
    }
    keep <- which(keep);
    statesSS <- states[keep];
    ySS <- y[keep];
    xSS <- x[keep];

    # Identify target loci with this state
    if (is.na(state)) {
      keep <- is.na(statesOut);
    } else {
      keep <- (statesOut == state);
    }
    keep <- which(keep);
    xOutSS <- xOut[keep];

    verbose && enter(verbose, "Kernel smoothing");
    verbose && cat(verbose, "Arguments:");
    args <- list(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, args);
    yOutSS <- kernelSmoothing(y=ySS, x=xSS, xOut=xOutSS, ...);
    verbose && str(verbose, yOutSS);
    verbose && exit(verbose);

    yOut[keep] <- yOutSS;
    verbose && exit(verbose);
  } # for (ss ...)
  verbose && str(verbose, yOut);

  verbose && enter(verbose, "Creating result object");
  # Allocate results of the correct size
  res <- newInstance(this, nrow=length(xOut));

  res$x <- xOut;
  res <- setSignals(res, yOut);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # kernelSmoothingByState()



setMethodS3("binnedSmoothingByState", "SegmentedGenomicSignalsInterface", function(this, from=xMin(this), to=xMax(this), by=NULL, length.out=NULL, byCount=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  x <- getPositions(this);
  # Argument 'from' & 'to':
  if (is.null(from)) {
    from <- min(x, na.rm=TRUE);
  } else {
    from <- Arguments$getInteger(from);
  }
  if (is.null(to)) {
    to <- max(x, na.rm=TRUE);
  } else {
    to <- Arguments$getInteger(to, range=c(from, Inf));
  }

  # Arguments 'by' & 'length.out':
  if (is.null(by) & is.null(length.out)) {
    throw("Either argument 'by' or 'length.out' needs to be given.");
  }
  if (!is.null(by)) {
    by <- Arguments$getNumeric(by, range=c(0,to-from));
  }
  if (!is.null(length.out)) {
    length.out <- Arguments$getInteger(length.out, range=c(1,Inf));
  }

  # Argument 'byCount':
  byCount <- Arguments$getLogical(byCount);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Binning data set");
  verbose && cat(verbose, "by: ", by);
  verbose && cat(verbose, "length.out: ", length.out);
  verbose && cat(verbose, "byCount: ", byCount);

  verbose && enter(verbose, "Find target positions");

  if (byCount) {
    verbose && enter(verbose, "By count");
    res <- sort(this); # sort() returns a clone():d object. /HB 2012-03-01
    resOut <- binnedSmoothing(res, by=by, length.out=length.out,
                              byCount=TRUE, verbose=less(verbose, 5));
    xOut <- resOut$x;
    # Not needed anymore
    resOut <- NULL;
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "By position");
    # Target 'x':
    if (!is.null(by)) {
      xOut <- seq(from=from, to=to, by=by);
    } else {
      xOut <- seq(from=from, to=to, length.out=length.out);
    }
    verbose && exit(verbose);
  } # if (byCount)

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);
  # Sanity check
  xOut <- Arguments$getNumerics(xOut);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating result set");
  # Allocate results of the correct size
  res <- newInstance(this, nrow=length(xOut));

  # Target 'x' and 'y':
  res$x <- xOut;
  ys <- rep(as.double(NA), times=length(xOut));
  res <- setSignals(res, ys);
  # Not needed anymore
  ys <- NULL;

  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying output states (may drop short regions)");
  statesOut <- getStates(res);
  verbose && cat(verbose, "statesOut:");
  verbose && str(verbose, statesOut);
  uniqueStates <- unique(statesOut);
  uniqueStates <- sort(uniqueStates, na.last=TRUE);
  verbose && cat(verbose, "Unique output states:");
  verbose && str(verbose, uniqueStates);
  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up source signals");
  if (byCount) {
    # Adding ordering along genome
    gs <- clone(this);
    gs <- sort(gs);
    gs$xOrder <- seq_len(nbrOfLoci(gs));
  } else {
    gs <- this;
  }
  verbose && print(verbose, gs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Binning (target) state by state
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ss in seq_along(uniqueStates)) {
    state <- uniqueStates[ss];
    verbose && enter(verbose, sprintf("State #%d ('%d') of %d",
                                        ss, state, length(uniqueStates)));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify target loci with this state
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Extracting subset of (target) loci with this signal state");
    idxsOut <- which(is.element(statesOut, state));
    resSS <- extractSubset(res, idxsOut, verbose=less(verbose,50));
    verbose && print(verbose, resSS);
    xOutSS <- getPositions(resSS);
    verbose && str(verbose, xOutSS);
    verbose && exit(verbose);
    # Nothing to do? [Should actually never happen!]
    nbrOfBins <- nbrOfLoci(resSS);
    # Not needed anymore
    resSS <- NULL;
    if (nbrOfBins == 0) {
      verbose && cat(verbose, "No bins. Skipping state.");
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identifying source loci with this state
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Extracting subset of (source) loci with this signal state");
    gsSS <- extractSubsetByState(gs, states=state, verbose=less(verbose,50));
#    verbose && print(verbose, gsSS);
    verbose && exit(verbose);
    # Nothing to do?
    if (nbrOfLoci(gsSS) == 0) {
      verbose && cat(verbose, "No extracted loci. Skipping state.");
      verbose && exit(verbose);
      next;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extending bins to equal count sizes?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (byCount) {
      verbose && enter(verbose, "Extending bins to have equal number of loci");
      verbose && cat(verbose, "xOrder:");
      verbose && str(verbose, gsSS$xOrder);
      nbrOfLociToAdd <- (gsSS$xOrder[1] %% by) - 1;
      verbose && cat(verbose, "Number of loci to add: ", nbrOfLociToAdd);
      gsSS$xOrder <- NULL;

      if (nbrOfLociToAdd > 0) {
        # Allocate the correct size
        nbrOfLoci2 <- nbrOfLoci(gsSS) + nbrOfLociToAdd;
        gsSS2 <- newInstance(gsSS, nrow=nbrOfLoci2);

        fields <- getDefaultLocusFields(gsSS, translate=TRUE);
        for (ff in fields) {
          values <- gsSS[[ff]];
          if (is.element(ff, "w")) {
            naValue <- 0;
          } else if (is.element(ff, "x")) {
            naValue <- 0;
          } else {
            naValue <- NA;
          }
          storage.mode(naValue) <- storage.mode(values);
          naValues <- rep(naValue, times=nbrOfLociToAdd);
          values <- c(naValues, values);
          gsSS2[[ff]] <- values;
        } # for (ff ...)
        gsSS <- gsSS2;
        # Not needed anymore
        gsSS2 <- NULL;
      } # if (nbrOfLociToAdd > 0)
      verbose && exit(verbose);
    } # if (byCount)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bin loci of this state towards target loci (of the same state)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Binned smoothing of temporary object");
    if (byCount) {
      xOutSSRank <- seq(from=(1+by)/2, by=by, length.out=nbrOfBins);
      verbose && cat(verbose, "Arguments:");
      args <- list(xOut=xOutSSRank, by=by, byCount=byCount, ...);
      verbose && str(verbose, args);
      resSS <- binnedSmoothing(gsSS, xOut=xOutSSRank, by=by, byCount=byCount, ...);
      resSS$x <- xOutSS;
    } else {
      verbose && cat(verbose, "Arguments:");
      args <- list(xOut=xOutSS, by=by, byCount=byCount, ...);
      verbose && str(verbose, args);
      resSS <- binnedSmoothing(gsSS, xOut=xOutSS, by=by, byCount=byCount, ...);
    } # if (byCount)
    # Not needed anymore
    gsSS <- args <- NULL;
    verbose && print(verbose, resSS);
    verbose && exit(verbose);

    # Sanity check
    stopifnot(length(idxsOut) == nbrOfLoci(resSS));

    ys <- getSignals(res);
    ys[idxsOut] <- getSignals(resSS);
    res <- setSignals(res, ys);
    # Not needed anymore
    resSS <- ys <- NULL;

    verbose && exit(verbose);
  } # for (ss ...)

  verbose && print(verbose, res);
  verbose && exit(verbose);

  res;
}) # binnedSmoothingByState()




setMethodS3("plot", "SegmentedGenomicSignalsInterface", function(x, ..., col=getStateColors(x)) {
  NextMethod("plot", col=col);
})

setMethodS3("points", "SegmentedGenomicSignalsInterface", function(x, ..., col=getStateColors(x)) {
  NextMethod("points", col=col);
})



############################################################################
# HISTORY:
# 2013-12-12
# o Added argument 'na.rm' to getUniqueStates().
# 2013-12-11
# o Added findChangePointsByState() for SegmentedGenomicSignalsInterface.
# 2012-03-01
# o Now binnedSmoothingByState() works also when RawGenomicSignals
#   extends data.frame.
# 2009-06-30
# o Added support for argument 'byCount' to binnedSmoothingByState() of
#   SegmentedCopyNumbers.  It is rather complex how it works, but we tried
#   to immitate how it works with byCount=FALSE.
# 2009-06-10
# o Added setStates().
# o Created SegmentedGenomicSignalsInterface class from SegmentedCopyNumbers.R.
# 2009-05-16
# o Now all methods of SegmentedRawCopyNumbers() coerce numerics only if
#   necessary, i.e. it keeps integers if integers, otherwise to doubles.
#   This is a general design of aroma.* that saves some memory.
# 2009-04-06
# o Now binnedSmoothingByState() of SegmentedCopyNumbers uses
#   extractSubsetByState() and then binnedSmoothing() on that object.
#   This makes the code slightly less redundant.
# 2009-02-19
# o Adopted to make use of new RawGenomicSignals.
# 2009-02-16
# o Now getStates() also passes the optional 'name' field to the "truth"
#   function.
# 2009-02-08
# o Added getUniqueStates().
# 2009-02-07
# o Added extractSubsetByState().
# o Created.
############################################################################
