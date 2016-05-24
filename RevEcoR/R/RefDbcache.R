#' Reference data for global metabolic construction 
#'  
#' The reference metabolic pathway data contains KOs, substrates and products, as well as 
#' a constructed reference global network, which used for metabolic network
#' reconstruction
#' 
#'     
#' @format The format is: List of 7 KO, substrate, product, user, date, version, 
#' reference network
#' 
#' @name RefDbcache
#' @rdname RefDbcache
#' 
#' @details Information this dataset is involved:
#' 
#' \itemize{
#'   \item KO, all KEGG orthlogy enties in KEGG metabolic pathways.
#'   \item substrate, substrate of enzymatic reactions in all KEGG metabolic pathways.
#'   \item product, product of enzymatic reactions in all KEGG metabolic pathways.
#'   \item user who download this data.
#'   \item date, the date this data is downloaded.
#'   \item version, R version used to obtained it.
#'   \item network, the global network which is reconstructed based on all the metabolites.
#' } 
#' 
#' @references \url{https://www.bioconductor.org/packages/release/bioc/html/mmnet.html}

NULL 