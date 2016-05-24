#' Human signal peptides cleavage sites
#' 
#' A set of 648 cleavage sites and 648 parts of mature proteins shortly after cleavage sites 
#' derived from human proteome.
#' 
#' @name human_cleave
#' @docType data
#' @details Each peptide in the data set is nine amino acid residues long. In case of cleavage sites,
#' the clevage is located between fifth and sixth peptide.
#' The non-cleavage sites are parts of mature proteins starting five positions after cleavage site.
#' @format A data frame with 1296 observations on the following 10 variables. Columns from
#' \code{P1} to \code{P9} describes positions in an extracted peptide. \code{tar} is a target vector. It
#' has value 1 if a peptide is a cleavage site and 0 if not.
#' @note Amino acid residues were recoded as integers.
#' @keywords datasets
#' @source \href{http://www.uniprot.org/}{UniProt}
#' @examples
#' 
#' data(human_cleave)
#' table(human_cleave[, 1])
#' 
NULL