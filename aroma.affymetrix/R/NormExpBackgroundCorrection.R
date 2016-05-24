###########################################################################/**
# @RdocClass NormExpBackgroundCorrection
#
# @title "The NormExpBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the normal exponential background correction model.
#  Estimators of the \pkg{limma} package is used.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "LimmaBackgroundCorrection".}
#   \item{method}{The estimator used, cf. argument \code{normexp.method}
#     of @see "limma::backgroundCorrect" in \pkg{limma} for more details.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   Internally, @see "limma::backgroundCorrect" is used.
# }
#*/###########################################################################
setConstructorS3("NormExpBackgroundCorrection", function(..., method=c("rma", "saddle", "mle")) {
  # Argument 'method':
  if (!is.null(method)) {
    method <- match.arg(method);
  }

  extend(LimmaBackgroundCorrection(...), "NormExpBackgroundCorrection",
    .method = method
  );
})



setMethodS3("getAsteriskTags", "NormExpBackgroundCorrection", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Drop added 'normexp' tag
  tags <- setdiff(tags, "normexp");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "NormExpBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Overload the 'args' for calling limma
  args <- params$args;
  args$method <- "normexp";
  args$normexp.method <- this$.method;
  params$args <- args;

  params;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2009-04-16
# o Added redundancy test that the default NormExpBackgroundCorrection
#   gives identical results to older RmaBackgroundCorrection.
# o Created NormExpBackgroundCorrection which extends more generic
#   LimmaBackgroundCorrection.
############################################################################
