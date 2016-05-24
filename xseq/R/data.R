#' TCGA AML SNP6.0 copy number alteration data
#'
#' A dataset containing part of the The Cancer Genome Atlas acute myeloid leukemia 
#' Affymetrix SNP6.0 array copy number alteration log2 ratios
#'
#' @format A matrix containing the copy number log2 ratios of 454 genes in 197 patients:
#' \itemize{
#'   \item{Row names are patient identifiers}
#'   \item{Colume names are official HGNC gene symbols}
#' }
#' @source \url{https://www.synapse.org/#!Synapse:syn300013}
"cna.logr"


#' TCGA AML SNP6.0 GISTIC copy number alteration calls
#'
#' A dataset containing part of the The Cancer Genome Atlas acute myeloid leukemia 
#' Affymetrix SNP6.0 array copy number alteration calls from GISTIC
#'
#' @format A matrix containing the GISTIC copy number calls of 454 genes in 197 patients:
#' \itemize{
#'   \item{Row names are patient identifiers}
#'   \item{Colume names are official HGNC gene symbols}
#' }
#' Each element of the matrix is coded:
#' \itemize{
#'   \item -2, homozygous deletions
#'   \item -1, hemizygous deletions
#'   \item 0, neutral
#'   \item 1, gain
#'   \item 2, amplifications
#'   }
#' @source \url{https://www.synapse.org/#!Synapse:syn300013}
"cna.call"


#' TCGA AML SNP6.0 gene expression data
#'
#' A dataset containing part of the The Cancer Genome Atlas acute myeloid leukemia RNA-seq gene expression data
#'
#' @format A matrix containing the expression of 454 genes in 197 patients:
#' \itemize{
#'   \item{Row names are patient identifiers}
#'   \item{Colume names are official HGNC gene symbols}
#' }
#' @source \url{https://www.synapse.org/#!Synapse:syn300013}
"expr"


#' TCGA AML somatic mutation data
#'
#' A dataset containing the The Cancer Genome Atlas acute myeloid leukemia somatic mutation data
#'
#' @format A data frame with 2311 rows and 12 variables:
#' \itemize{
#'   \item sample. character, patient identifier
#'   \item hgnc_symbol. character, official HGNC gene symbols
#'   \item variant_type. character, mutation type, can be either "HOMD", "HLAMP", 
#'   "MISSENSE", "NONSENSE", "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
#'   "STARTGAINED", "SYNONYMOUS", "OTHER", "FUSION", "COMPLEX"
#'   \item ...
#' }
#' @source \url{https://www.synapse.org/#!Synapse:syn1729383}
"mut"


#' A networks containing gene associations
#'
#' A list of gene interactions
#'
#' @format A list with 16 elements:
#' \itemize{
#'   \item{List names are official HGNC gene symbols}
#'   \item{Each element of the list is a vector, and the name of the vector is an official HGNC gene symbol. 
#'   Each element of the vector is a real number between 0 and 1, 
#'   representing the assocation strongth between two genes. 
#'   The names of the elements of the vector are also official HGNC gene symbols.}
#' }
"net"

