#' Soybean genealogical data
#' 
#' This data contains information on copy number variants, single nucleotide polymorphisms, protein content, and yield, of soybeans. The available data consists of a data frame structure that contains 412 direct child-parent relationships between pairs of soybean varieties. These data were collected from field trials, genetic studies, and United States Department of Agriculture (USDA) bulletins, and date as early as the first decade of the 1900s.

#' @name sbGeneal
#' @title Soybean data
#' @description This data set contains soy bean genealogical information maintained by the United States Department of Agriculture to be used by plant breeders, geneticists, bioinformaticians, pathologists,and many other research workers.
#' @docType data
#' @format a \code{RData} instance, 1 row per each child-parent relationship between soybean varieties
#' @details \itemize{
  #' \item child name of child soybean variety
  #' \item year child variety was introduced
  #' \item yield protein yield
  #' \item year.imputed whether or not the introduced year of the child variety was imputed
  #' \item parent name of parent soybean variety
  #' }
  #'
#' @docType data
#' @keywords datasets
#' @name sbGeneal
#' @usage data(sbGeneal)
#' @format A data frame with 412 rows and 7 variables
#' @references
#' "Pedigrees of soybean cultivars released in the United States and Canada." Theodore Hyivitz, C.A. Newell, S.G. Carmer. College of Agriculture, University of Illinois at Urbana-Champaign (1977).
NULL
