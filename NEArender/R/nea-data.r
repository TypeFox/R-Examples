#' NEArender
#' 
#' TCGA point mutations from the glioblastoma multiforme cohort
#'
#' @docType data
#' @name tcga.gbm
#' @usage data(tcga.gbm)
#' @keywords datasets
#' @format A matrix  with genes as rows (9660 genes with mutations reported in at least one sample) and TCGA samples as columns (291 samples). Wild type (=reference) states are denoted with NA.
#' @examples
#' table(is.na(tcga.gbm))
NULL


#' Example functional gene sets (FGS) from GO and KEGG
#'
#' @docType data
#' @name can.sig.go
#' @usage data(can.sig.go)
#' @keywords datasets
#' @format An object, i.e. a list with FGS IDs as names and vectors of gene IDs as member (34 FGSs in total). See see \code{\link{import.gs}}
#'
#' @examples
#' head(can.sig.go)
#'
NULL

#' A Fantom5 transcriptomes in 43 carincinoma cell samples
#'
#' @docType data
#' @name fantom5.43samples
#' @keywords datasets
#' @usage data(fantom5.43samples)
#' @format An expression data matrix with the genes as rows and columns as cell line samples
#'
#' @examples
#' colnames(fantom5.43samples)
#'
NULL

#' A more compact alternative global network for the network enrichment analysis
#'
#' @docType data
#' @name net.kegg
#' @keywords datasets
#' @usage data(net.kegg)
#' @format A list, see \code{\link{import.net}}
#'
#' @examples
#' head(net.kegg)
#'
"net.kegg"
