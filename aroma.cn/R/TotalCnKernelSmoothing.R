###########################################################################/**
# @RdocClass TotalCnKernelSmoothing
#
# @title "The TotalCnKernelSmoothing class"
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
#  \item{kernel}{A @character string specifying the type of kernel
#     to be used.}
#  \item{bandwidth}{A @double specifying the bandwidth of the smoothing.}
#  \item{censorH}{A positive @double specifying the bandwidth threshold
#     where values outside are ignored (zero weight).}
#  \item{robust}{If @TRUE, a robust smoother is used, otherwise not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("TotalCnKernelSmoothing", function( ..., kernel=c("gaussian", "uniform"), bandwidth=50e3, censorH=3, robust=FALSE) {
  # Argument 'kernel':
  kernel <- match.arg(kernel);

  # Argument 'bandwidth':
  bandwidth <- Arguments$getDouble(bandwidth, range=c(0,Inf));

  # Arguments 'censorH':
  censorH <- Arguments$getDouble(censorH, range=c(0,Inf));

  # Arguments 'robust':
  robust <- Arguments$getLogical(robust);

  extend(TotalCnSmoothing(...), "TotalCnKernelSmoothing",
    .kernel = kernel,
    .bandwidth = bandwidth,
    .censorH = censorH,
    .robust = robust
  );
})


setMethodS3("getParameters", "TotalCnKernelSmoothing", function(this, ...) {
  params <- NextMethod("getParameters");
  params$kernel <- this$.kernel;
  params$bandwidth <- this$.bandwidth;
  params$censorH <- this$.censorH;
  params$robust <- this$.robust;
  params;
}, protected=TRUE);


setMethodS3("getAsteriskTags", "TotalCnKernelSmoothing", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class-specific tags

  params <- getParameters(this);

  # We put the bandwidth tag before the kernel one for
  # backward compatibility reason. /HB 2011-12-15
  # Parameter 'bandwidth'
  bandwidthTag <- sprintf("H=%.1fkb", params$bandwidth/1e3);
  tags <- c(tags, bandwidthTag);

  # "Parameter" 'kernel'
  kernelTag <- params$kernel;
  tags <- c(tags, kernelTag);

  # Parameter 'robust'
  if (params$robust)
    tags <- c(tags, "robust");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("smoothRawCopyNumbers", "TotalCnKernelSmoothing", function(this, rawCNs, target, ..., verbose=FALSE) {
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
  args <- c(list(xOut=target$xOut), params, h=params$bandwidth, ...);

  # Keep only known arguments
  knownArguments <- names(formals(colKernelSmoothing.matrix));
  keep <- is.element(names(args), knownArguments);
  args <- args[keep];

  args <- c(list(rawCNs), args);

  verbose && cat(verbose, "Calling kernelSmoothing() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 20);
  smoothCNs <- do.call("kernelSmoothing", args=args);

  verbose && exit(verbose);

  smoothCNs;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2011-12-15
# o Moved argument 'bandwidth' to TotalCnKernelSmoothing.
# 2009-02-08
# o Created from TotalCnSmoothing.R.
# 2009-01-26
# o Adopted to the new AromaUnitTotalCnBinarySet.
# o Added Rdoc comments.
# 2008-05-23
# o Created.
############################################################################
