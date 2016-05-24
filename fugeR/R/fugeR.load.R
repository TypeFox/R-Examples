######################################################################################
# fugeR.load
#
#'   Load a fuzzy system.
#'   
#'   Load a fuzzy system saved into a file with \code{fugeR.save}
#'   
#'   @param file        [\"\"] A character string naming a file.
#'
#'   @examples
#'   ##
#'   ##
#'   \dontrun{
#'      fis <- fugeR.run (
#'                  In,
#'                  Out,
#'                  generation=100,
#'                  population=200,
#'                  elitism=40,
#'                  verbose=TRUE,
#'                  threshold=0.5,
#'                  sensiW=1.0,
#'                  speciW=1.0,
#'                  accuW=0.0,
#'                  rmseW=1.0,
#'                  maxRules=10,
#'                  maxVarPerRule=2,
#'                  labelsMf=2
#'              )
#'      
#'      fugeR.save( fis, file=\'./myFis.R\' )
#'      
#'      savedFis <- fugeR.load( file=\'./myFis.R\' )
#'   }
#'
#' @seealso  \code{\link{fugeR.save}}
#' 
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.load <-
function(file="") {
    #check if a path was specified
    if(file == "") {
        stop("File name can't be empty")
    }
    dget(file)
} 
