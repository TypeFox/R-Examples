#' Creating an updated ecogen object by removing results of the slot OUT
#' 
#' @param eco Object of class "ecogen".
#' @param ... Objects to remove from eco, typed without quotations. 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' variog <- eco.variogram(eco[["P"]][, 1], eco[["XY"]])
#' 
#' # Assignation of values can be made with the corresponding accessors,
#' # using the generic notation of EcoGenetics 
#' # (<ecoslot.> + <name of the slot> + <name of the object>).
#' # See help("EcoGenetics accessors")
#' 
#' ecoslot.OUT(eco) <- variog         
#' we.are.numbers <- c(1:10)
#' we.are.characters <- c("John Coltrane", "Charlie Parker")
#' ecoslot.OUT(eco) <- list(we.are.numbers, we.are.characters)
#' ecoslot.OUT(eco)
#' eco <- eco.remove(eco, we.are.numbers)
#' ecoslot.OUT(eco)
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export


setGeneric("eco.remove", 
           
           function(eco, ...) {
             
             res.names <- as.character(match.call())
             res.names <- res.names[-c(1:2)]
             del <- (names(eco@OUT)) %in% res.names
             eco@OUT <- eco@OUT[!del]
             eco
           })
