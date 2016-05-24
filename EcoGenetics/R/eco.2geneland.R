#' Creating input data for Geneland with an ecogen object
#' 
#' @description This function creates four data frames in the working 
#' directory (XY.txt, NAMES.txt, P.txt, G.txt) which can be loaded 
#' in Geneland.
#' @param eco Object of class "ecogen".
#' @param ncod Number of digits coding each allele
#'  (e.g., 1: x, 2: xx, 3: xxx, etc.).
#' @param ploidy Ploidy of the data.
#' @return XY.txt Matrix with coordinates.
#' @return NAMES.txt Matrix with row names.
#' @return P.txt Matrix with phenotypic data.
#' @return G.txt Matrix with genotypic data.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco.2geneland(eco, 1)
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


setGeneric("eco.2geneland", 
           function(eco, ncod = NULL, ploidy = 2) {
             
             
             write.table(eco@XY, "XY.txt", quote = FALSE,
                         row.names = FALSE, col.names = FALSE)
             
             write.table(rownames(eco@XY), "NAMES.txt", quote = FALSE, 
                         row.names = FALSE, col.names = FALSE)
             
             write.table(eco@P,"P.txt", quote = FALSE, row.names = FALSE, 
                         col.names = FALSE)
             
             write.table(int.loc2al(eco@G,  ncod = ncod,  ploidy = ploidy), "G.txt",
                         quote = FALSE, row.names = FALSE, col.names = FALSE)
             return("done!")
             
           })
