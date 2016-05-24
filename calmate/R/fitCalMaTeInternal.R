###########################################################################/**
# @RdocDocumentation fitCalMaTeInternal
# @alias fitCalMaTeV1
# @alias fitCalMaTeV2
# @alias fitCalMaTeMedians
#
# @title "Algorithms to fit the CalMaTe model for a single SNP"
#
# \description{
#  @get "title".
#  \emph{Note: These are internal functions of the package.
#        They should not be used elsewhere.}
# }
#
# \usage{
#   fitCalMaTeV1(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...)
#   fitCalMaTeV2(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...)
#   fitCalMaTeMedians(dataT, references, fB1=1/3, fB2=2/3,...)
# }
#
# \arguments{
#  \item{dataT}{A 2xI @numeric @matrix of allele specific copy numbers (ASCNs),
#     where 2 is the number alleles and I is the number of samples.}
#  \item{references}{A @integer @vector with elements in [1,I] specifying
#     which samples should be used as the reference set.}
#  \item{fB1, fB2}{Thresholds for calling genotypes AA, AB, BB from the
#     allele B fractions.}
#  \item{maxIter}{The maximum number of iterations without converging
#     before the algorithm quits.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a 2xI @numeric @matrix of calibrated ASCNs.
# }
#
# \section{Flavor v1}{
#   This is an early version (June 2010-January 2012) of the algorithm
#   described in [1].
# }
#
# \section{Flavor v2}{
#   This is the model and algorithm described in [1].
#
#   This version was introduced to decrease the number of
#   "artificial outliers" introduced by CalMaTe for some SNPs 
#   due to non-converging or wreak-havoc estimates of the SNP effects.
#   Flavor v2 differ from Flavor v1 as follows:
#   \itemize{
#    \item The estimation of the model parameters are now done solely
#          based on reference samples.  In previous versions, some of the
#          initial estimation steps were using also non-reference samples.
#    \item For a small number of SNPs, the main CalMaTe scheme for estimating
#          parameters would not converge or converge poorly.  For such SNPs
#          CalMaTe now falls back to using a plain median estimator,
#          i.e. \code{fitCalMaTeMedians()}.
#    \item The above fallback estimator is also used in cases where
#          all samples are indentified to be homozygous.
#   }
#
#   The \code{fitCalMaTeMedians()} method is used as a fallback method
#   by \code{fitCalMaTeV2()}.  It fits CalMaTe without using the
#   @see "MASS::rlm" function.
# }
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2012.Rd" \cr 
# }
#
# \seealso{
#   These functions are called by @see "calmateByThetaAB".
# }
#
# @keyword internal
#*/###########################################################################


###########################################################################
# HISTORY:
# 2012-02-21 [HB]
# o Added fitCalMaTeMedians() to the documentation.
# o Added documentation for Flavor v2.
# 2012-02-19 [HB]
# o This file holds the help for fitCalMaTeV1() and fitCalMaTeV2().
# o Created.
###########################################################################
