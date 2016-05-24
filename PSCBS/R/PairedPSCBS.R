###########################################################################/**
# @RdocClass PairedPSCBS
#
# @title "The PairedPSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PairedPSCBS is an object containing the results from the
#  Paired PSCBS method.
# }
#
# \usage{PairedPSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the Paired PSCBS results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   The @see "segmentByPairedPSCBS" method returns an object of this class.
# }
#*/###########################################################################
setConstructorS3("PairedPSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(PSCBS(fit=fit, ...), "PairedPSCBS");
})

setMethodS3("getLocusData", "PairedPSCBS", function(fit, ..., fields=c("asis", "full")) {
  # Argument 'fields':
  fields <- match.arg(fields);


  data <- NextMethod("getLocusData", fields="asis");


  if (fields == "full") {
    names <- colnames(data);

    # Genotype calls
    if (!is.element("muN", names)) {
      callNaiveGenotypes <- .use("callNaiveGenotypes", package="aroma.light");
      data$muN <- callNaiveGenotypes(data$betaN);
    }
    data$isHet <- (data$muN == 1/2);
    # BACKWARD COMPATIBILITY: If 'rho' does not exists, calculate
    # it on the fly from 'betaT'.
    # NOTE: This should give an error in the future. /HB 2013-10-25
    if (is.null(data$rho)) {
      data$rho <- 2*abs(data$betaT-1/2);
      data$rho[!data$isHet] <- NA_real_;
      warning("Locus-level DH signals ('rho') did not exist and were calculated from tumor BAFs ('betaT')");
    }
    data$c1 <- 1/2*(1-data$rho)*data$CT;
    data$c2 <- data$CT - data$c1;

    # TumorBoost BAFs
    if (!is.element("betaTN", names)) {
      normalizeTumorBoost <- .use("normalizeTumorBoost", package="aroma.light");
      data$betaTN <- normalizeTumorBoost(betaN=data$betaN, betaT=data$betaT, muN=data$muN);
    }
    data$rhoN <- 2*abs(data$betaTN-1/2);
    data$rhoN[!data$isHet] <- NA_real_;
    data$c1N <- 1/2*(1-data$rhoN)*data$CT;
    data$c2N <- data$CT - data$c1N;

    data$isSNP <- (!is.na(data$betaT) | !is.na(data$betaN));
    data$type <- ifelse(data$isSNP, "SNP", "non-polymorphic locus");

    # Labels
    data$muNx <- c("AA", "AB", "BB")[2*data$muN + 1L];
    data$isHetx <- c("AA|BB", "AB")[data$isHet + 1L];
  }

  data;
}, protected=TRUE) # getLocusData()



setMethodS3("resegment", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Resegmenting a ", class(fit)[1], " object");

  # Use the locus-level data of the PairedPSCBS object
  data <- getLocusData(fit);
  class(data) <- "data.frame";
  drop <- c("rho", "betaTN", "index");
  keep <- !is.element(colnames(data), drop);
  data <- data[,keep];
  verbose && str(verbose, data);

  verbose && cat(verbose, "Number of loci: ", nrow(data));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup arguments to be passed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Overriding default arguments");
  segFcnName <- "segmentByPairedPSCBS";
  segFcn <- getMethodS3(segFcnName, "default");

  # (a) The default arguments
  formals <- formals(segFcn);

  formals <- formals[!sapply(formals, FUN=is.language)];
  formals <- formals[!sapply(formals, FUN=is.name)];
  drop <- c("chromosome", "x", "w", "CT", "thetaT", "thetaN", "betaT", "betaN", "muN", "rho", "...");
  keep <- !is.element(names(formals), drop);
  formals <- formals[keep];

  # (b) The arguments used in previous fit
  params <- fit$params;
  keep <- is.element(names(params), names(formals));
  params <- params[keep];
  # Don't trust 'tbn'!  TODO. /HB 20111117
  params$tbn <- NULL;

  # (c) The arguments in '...'
  userArgs <- list(..., verbose=verbose);

  # (d) Merge
  args <- formals;
  args2 <- append(params, userArgs);
  for (kk in seq(along=args2)) {
    value <- args2[[kk]];
    if (!is.null(value)) {
      key <- names(args2)[kk];
      if (!is.null(key)) {
        args[[key]] <- value;
      } else {
        args <- append(args, list(value));
      }
    }
  } # for (key ...)
  verbose && str(verbose, args[names(args) != "verbose"]);

  verbose && enter(verbose, sprintf("Calling %s()", segFcnName));
  args <- append(list(data), args);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args[names(args) != "verbose"]);
  verbose && exit(verbose);

  fit <- do.call(segFcnName, args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # resegment()


setMethodS3("adjustPloidyScale", "PairedPSCBS", function(fit, scale, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Update locus-level data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  names <- c("CT");
  for (ff in names) {
    data[[ff]] <- scale * data[[ff]];
  }
  fit$data <- data;  ##  fit <- setLocusData(fit, data);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) Update segment-level data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit);

  # Adjust segment levels
  names <- grep("^(tcn|c1|c2)(Mean|_.*%)$", names(segs), value=TRUE);
  for (ff in names) {
    segs[[ff]] <- scale * segs[[ff]];
  }

  # Clear segment calls
  names <- c("lohCall", "ntcnCall");
  for (ff in names) {
    segs[[ff]] <- NULL;
  }
  fit$output <- segs; ## fit <- setSegments(fit, sets);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (c) Update parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- fit$params;
  fields <- c("copyNeutralStats", "deltaCN", "ntcnRange", "deltaLowC1");
  params[fields] <- NULL;
  fit$params <- params;

  fit;
}, protected=TRUE) # adjustPloidyScale()


##############################################################################
# HISTORY
# 2014-03-30
# o Update resegment() for PairedPSCBS to handle 'thetaT' and 'thetaN'.
# 2013-10-25
# o BUG FIX: The 'rho' signals returned by getLocusData(..., fields="full")
#   for PairedPSCBS would have values also for homozygote SNPs.
# 2013-03-08
# o Added getLocusData() for PairedPSCBS.
# 2012-04-21
# o CLEANUP: Moved getSegmentSizes() from PairedPSCBS to PSCBS.
# 2011-11-21
# o BUG FIX: resegment() was trying to call segmentByCBS() instead
#   of segmentByPairedPSCBS().
# 2011-11-17
# o Added resegment() for PairedPSCBS for easy resegmentation.
# 2011-10-02
# o CLEANUP: Moved print() and as.data.frame() to PSCBS.
# o Added Rdoc help.
# o Now the constructor of PairedPSCBS calls that of PSCBS.
# 2011-06-28
# o DOCUMENTATION: Added Rd help for as.data.frame() of PairedPSCBS.
# 2011-04-08
# o Added formal constructor for the PairedPSCBS class.
# o Created.
##############################################################################
