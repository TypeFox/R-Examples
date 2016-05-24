###########################################################################/**
# @RdocClass TotalCnBinnedSmoothing
#
# @title "The TotalCnBinnedSmoothing class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "TotalCnSmoothing".}
#  \item{robust}{If @TRUE, a robust smoother is used, otherwise not.}
# }
#
# \details{
#  Note that \code{dsS <- TotalCnBinnedSmoothing(ds, targetUgp=ugp)} where
#  \code{ugp <- getAromaUgpFile(ds)} returns a data set with an identical
#  set of loci as the input data set and identical signals as the
#  input ones, \emph{except} for loci with duplicated positions.  If all
#  loci have unique positions, the the output is identical to the input.
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("TotalCnBinnedSmoothing", function( ..., robust=FALSE) {
  # Arguments 'robust':
  robust <- Arguments$getLogical(robust);

  extend(TotalCnSmoothing(...), "TotalCnBinnedSmoothing",
    .robust = robust
  );
})


setMethodS3("getParameters", "TotalCnBinnedSmoothing", function(this, ...) {
  params <- NextMethod("getParameters");
  params$robust <- this$.robust;
  params;
}, protected=TRUE);


setMethodS3("getAsteriskTags", "TotalCnBinnedSmoothing", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class-specific tags

  params <- getParameters(this);

  # Parameter 'robust'
  if (params$robust)
    tags <- c(tags, "robust");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("smoothRawCopyNumbers", "TotalCnBinnedSmoothing", function(this, rawCNs, target, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing one set of copy numbers");
  verbose && print(verbose, rawCNs);

  # Setting up arguments
  params <- getParameters(this);
  params$FUN <- ifelse(params$robust, "median", "mean");
  params$robust <- NULL;
  args <- c(list(xOut=target$xOut), params, ...);

  # Keep only known arguments
  knownArguments <- names(formals(colBinnedSmoothing.matrix));
  keep <- is.element(names(args), knownArguments);
  args <- args[keep];

  args <- c(list(rawCNs), args);

  verbose && cat(verbose, "Calling binnedSmoothing() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 20);
  smoothCNs <- do.call("binnedSmoothing", args=args);

  verbose && exit(verbose);

  smoothCNs;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-08-26
# o DOCUMENTATION: Added a help section on binning with target loci being
#   identical to the input loci.
# 2012-01-16
# o Created from TotalCnKernelSmoothing.R.
############################################################################
