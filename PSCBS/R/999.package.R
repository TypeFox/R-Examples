#########################################################################/**
# @RdocPackage "PSCBS"
# @eval "rd <- gsub('\\alias{PSCBS}', '', rd, fixed=TRUE); ''"
#
# \description{
#   @eval "packageDescription('PSCBS')$Description".
#
#   This package should be considered to be in an alpha or beta phase.
#   You should expect the API to be changing over time.
# }
#
# \section{Installation and updates}{
#   To install this package, use \code{install.packages("PSCBS")}.
# }
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item Vignettes '\href{../doc/index.html}{Parent-specific copy-number segmentation using Paired PSCBS}' and '\href{../doc/index.html}{Total copy-number segmentation using CBS}'.
#     \item @see "segmentByCBS" - segments total copy-numbers, or any
#           other unimodal genomic signals, using the CBS method [3,4].
#     \item @see "segmentByPairedPSCBS" - segments allele-specific
#           tumor signal from a tumor \emph{with} a matched normal
#           using the Paired PSCBS method [1,2].
#     \item @see "segmentByNonPairedPSCBS" - segments allele-specific
#           tumor signal from a tumor \emph{without} a matched normal
#           using the Non-Paired PSCBS method adopted from [1,2].
#   }
# }
#
# \section{How to cite}{
#   Please use [1] and [2] to cite when using Paired PSCBS,
#   and [3] and [4] when using CBS.
#   When using Non-Paired PSCBS, please cite [1] and [2] as well.
# }
#
# @author
#
# \section{License}{
#  @eval "packageDescription('PSCBS')$License".
# }
#
# \references{
#  [1] @include "../incl/OlshenA_etal_2011.Rd" \cr
#  [2] @include "../incl/BengtssonH_etal_2010.Rd" \cr
#  [3] @include "../incl/OlshenVenkatraman_2004.Rd" \cr
#  [4] @include "../incl/VenkatramanOlshen_2007.Rd" \cr
# }
#*/#########################################################################

