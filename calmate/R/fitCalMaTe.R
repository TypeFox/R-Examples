###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTe
# @alias fitCalMaTe
#
# @title "Calibrates SNP loci according to the CalMaTe method"
#
# \description{
#  @get "title".
#  \emph{Note: This is an internal function of the package, which is kept
#   only kept to provide easy access to the internal fit functions.
#   It it actually not elsewhere in the package, and should nor by others.}
# }
#
# @synopsis
#
# \arguments{
#  \item{dataT}{A 2xI @numeric @matrix of allele specific copy numbers (ASCNs),
#     where 2 is the number alleles and I is the number of samples.}
#  \item{references}{A @integer @vector with elements in [1,I] specifying
#     which samples should be used as the reference set.}
#  \item{...}{Additional arguments passed to the internal fit functions.}
#  \item{flavor}{A @character string specifying which internal fit function
#     (flavor of the CalMaTe algorithm) to use for fitting the model.}
# }
#
# \value{
#   Returns a 2xI @numeric @matrix of calibrated ASCNs.
# }
#
# \section{Flavors}{
#   For backward compatibility, we try to keep all major versions of
#   the CalMaTe algorithm available.  Older versions can be used by
#   specifying argument \code{flavor}.
#   For more information about the different flavors, 
#   see @see "fitCalMaTeInternal".
# }
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2012.Rd" \cr 
# }
#
# \seealso{
#   For further information on the internal fit functions, 
#   see @see "fitCalMaTeInternal".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitCalMaTe", "matrix", function(dataT, references, flavor=c("v2", "v1"), ...) {
  if (flavor == "v2") {
    res <- fitCalMaTeV2(dataT, references=references, ...);
  } else if (flavor == "v1") {
    res <- fitCalMaTeV1(dataT, references=references, ...);
  } else {
    throw("Unknown algorithm flavor: ", flavor);
  }
  res;
}, protected=TRUE) # fitCalMaTe()


###########################################################################
# HISTORY:
# 2012-02-19 [HB]
# o Note that fitCalMaTe() for matrices is only kept for an easy access
#   to internal SNP CalMaTe estimators.  It is not used by 
#   calmateByThetaAB() nor calmateByTotalAndFracB(), which call the
#   internal fit functions directly.
# o BACKWARD COMPATIBILITY: Added argument 'flavor' to fitCalMaTe()
#   to be able to use previous versions of CalMaTe model estimators.
#   Flavor "v2" was introduced 2011-12-05.
#
# For the history of the CalMaTe algorithm, see fitCalMaTeV*.R.
###########################################################################
