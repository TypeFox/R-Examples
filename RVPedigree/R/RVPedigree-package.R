#' RVPedigree: A package for region-based genetic
#' association tests.
#'
#' The RVPedigree package contains three methods for
#' doing region-based association tests. It contains methods to 
#' perform autosomal rare variant association analyses with a normally-distributed or non 
#' normally distributed quantitative phenotype, 
#' specifically for the situation when there are related individuals or families in the data. 
#' The methods included are:
#' \itemize{
#'   \item ASKAT (exact test)
#'   \item Normalized ASKAT (exact test)
#'   \item VC-C1 (exact test)
#'   \item VC-C2 (permutation test of VC-C1 statistic)
#'   \item VC-C3 (permutation test of VC-C1 statistic)
#' }
#'
#'
#' @section RVPedigree functions:
#' The following functions are the user-visible functions of the
#'     package:
#' \itemize{
#' \item \code{\link{RVPedigree}}: main function that runs the
#'     selected test (ASKAT by default) genome-wide. If you only want
#'     to run a given test on a single genomic region, use one of the
#'     following functions.
#' \item \code{\link{ASKAT.region}}: Runs the ASKAT test on a given
#'     genomic region.
#' \item \code{\link{NormalizedASKAT.region}}: Runs the Normalized
#'     ASKAT test on a given genomic region.
#' \item \code{\link{VCC1.region}}: Runs the VC-C1 test on a given
#'     genomic region.
#' \item \code{\link{VCC2.region}}: Runs the VC-C2 test on a given
#'     genomic region.
#' \item \code{\link{VCC3.region}}: Runs the VC-C3 test on a given
#'     genomic region.
#' \item \code{\link{Estim.H0.ASKAT}}: Estimates the null model for
#'     the ASKAT test. The result of this function can be passed to
#'     the corresponding \code{\link{ASKAT.region}} function to save
#'     computation time in case multiple genomics regions are to be
#'     analyzed.
#' \item \code{\link{Estim.H0.NormalizedASKAT}}: Estimates the null
#'     model for the Normalized ASKAT test. The result of this
#'     function can be passed to the corresponding
#'     \code{\link{NormalizedASKAT.region}} function to save
#'     computation time in case multiple genomics regions are to be
#'     analyzed.
#' \item \code{\link{Estim.H0.VCC}}: Estimates the null model for the
#'     VC-C1, VC-C2 and VC-C3 tests. The result of this function can
#'     be passed to the corresponding \code{VCC?.region} function to
#'     save computation time in case multiple genomics
#'     regions are to be analyzed.
#'
#' \item \code{\link{GetRelMatrix}}: calculates the relationship matrix
#'     (twice the kinship matrix) based on various types of input
#'     data.
#'
#' \item \code{\link{readMapFile}}: reads a genetic map file (e.g. from
#'     Plink data) and creates the correct data frame to pass on to
#'     the various \code{*.region} functions.
#'
#' \item \code{\link{Normality.test}}: function to test for the
#'     normality of the phenotype data.
#' }
#'
#' @docType package
#' @name RVPedigree-package
#'
#' @seealso \code{GenABEL}, \code{snpStats}
#' @keywords package
NULL
#> NULL
