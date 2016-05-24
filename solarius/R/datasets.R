#' dat30 data set adapted from multic R package
#'
#' 29 first families were selected from the complete data set of 12000 individuals.
#' For a resulted subset of 174 individuals,
#' a hundred of synthetic SNPs were randomly generated.
#' Annotation information also was generated, 
#' mainly in order to plot the association results with Manhattan plot.
#' 
#' Two simulated phenotypes possess a high genetic correlation.
#'
#' @format 
#' (Phenotypes) A data frame \code{dat30} with 174 rows and 10 variables:
#' \describe{
#'   \item{famid}{Family ID (29 unique ids).}
#'   \item{id}{Individual ID.}
#'   \item{fa}{Father ID.}
#'   \item{mo}{Mother ID.}
#'   \item{sex}{Individual gender (1 - male, 2 - female).}
#'   \item{affect}{Affected status (1 - unaffected, 2 - affected).}
#'   \item{class}{Class label.}
#'   \item{trait1}{Simulated phenotype 1.}
#'   \item{trait2}{Simulated phenotype 2.}
#'   \item{age}{Age.}
#' }
#'
#' (Genotypes as covariates) A matrix \code{genocovdat30} with 174 rows and 100 columns.
#' Row names are IDs of individuals, column names are names of SNPs.
#'
#' (Annotation) A data frame \code{mapdat30} with 100 rows and 4 variables:
#' \describe{
#'   \item{SNP}{SNP name.}
#'   \item{chr}{Chromosome.}
#'   \item{pos}{Position in bp.}
#'   \item{gene}{Gene.}
#' }
#'
#' @source \url{https://cran.r-project.org/package=multic}
#'
#' @name dat30
#' @rdname dat30
#'
#' @examples
#' data(dat30)
#'
#' str(dat30)
#'
#' plotPed(dat30, 2) # plot the pedigree tree for family #2
#'
#' \dontrun{
#' kin2 <- solarKinship2(dat30)
#' plotKinship2(kin2) 
#' plotKinship2(kin2[1:30, 1:30])
#'
#' }
#'
NULL


#' @name genocovdat30
#' @rdname dat30
#'
#' @examples
#'
#' str(genocovdat30)
#'
#' genocovdat30[1:5, 1:5]
#'
NULL

#' @name mapdat30
#' @rdname dat30
#'
#' @examples
#'
#' str(mapdat30)
#'
#' head(mapdat30)
#'
NULL


#' dat50 data set adapted from FFBSKAT R package
#'
#' A mixture of unrelated and related individuals 
#' was originally simulated in FFBSKAT R package to test methods of the variant-collapsing approach.
#' 50 synthetic SNPs were generated for the association study.
#' A custom kinship matrix is used to express the relationships among individuals.
#' This data set is used here to test the ability of SOLAR 
#' to work with a custom kinship matrix in both polygenic and association analyses.
#'
#'
#' The genotypes are coded in the format such as 1/1, 1/2 and 2/2.
#'
#' In addition to the original data set from FFBSKAT R package,
#' a matrix of covariates was derived from the genotype data according to the additive model.
#'
#' @format (Phenotypes) A data frame \code{phenodata} with 66 rows and 4 variables:
#' \describe{
#'   \item{id}{Individual ID.}
#'   \item{sex}{Individual gender (0 - male, 1 - female).}
#'   \item{age}{Age.}
#'   \item{trait}{Simulated phenotype.}
#' }
#'
#' (Kinship) A square matrix \code{kin} with 66 rows and 66 columns.
#'
#' (Genotypes) A matrix \code{genodata} with 66 rows and 50 columns.
#'
#' (Genotypes as covariates) A matrix \code{genocovdata} with 66 rows and 50 columns.
#'
#' (Annotation) A data frame \code{snpdata} with 100 rows and 4 variables:
#' \describe{
#'   \item{name}{SNP name.}
#'   \item{chrom}{Chromosome.}
#'   \item{position}{Position in bp.}
#'   \item{gene}{Gene.}
#' }
#'
#' @name phenodata
#' @rdname dat50
#'
#' @source \url{http://mga.bionet.nsc.ru/soft/FFBSKAT/}
#'
#' @examples
#' data(dat50)
#'
#' str(phenodata)
#'
NULL

#' @name kin
#' @rdname dat50
#'
#' @examples
#'
#' plotKinship2(2*kin)
#'
NULL

#' @name genodata
#' @rdname dat50
#'
#' @examples
#'
#' str(genodata)
#'
#' genodata[1:5, 1:5]
#'
NULL


#' @name genocovdata
#' @rdname dat50
#'
#' @examples
#'
#' str(genocovdata)
#'
#' genocovdata[1:5, 1:5]
#'  
#' # compare with the genotypes
#' genodata[1:5, 1:5]
#'
NULL


#' @name snpdata
#' @rdname dat50
#'
#' @examples
#'
#' str(snpdata)
#'
#' head(snpdata)
#'
NULL


