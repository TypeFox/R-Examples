######################################################################################
# fugeR.save
#
#'   Save a fuzzy system into a file.
#'   
#'   @param fuzzySystem [NULL] The fuzzy system to save.
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
#' @seealso  \code{\link{fugeR.load}}
#' 
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.save <-
function(fuzzySystem=NULL, file="") {
    #check if it's really a fuzzy system
    if(class(fuzzySystem) != "fuzzySystem") {
        stop("Wrong object, fuzzySystem is not a valid fuzzy System")
    }
    #check if a path was specified
    if(file == "") {
        stop("File name can't be empty")
    }
    dput(fuzzySystem, file)
}