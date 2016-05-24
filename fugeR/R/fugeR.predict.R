######################################################################################
# fugeR.predict
# TO DO: Change to make the prediction even if the dataset has not the same number of column.
#
#'   Compute the prediction of a fuzzy system for the given input data.
#'   
#'   @param fuzzySystem [NULL] The fuzzy system to use for computing the prediction.
#'   @param dataset     [NULL] The data to use.
#'
#'   @return  \code{prediction}, A data.frame containing the predictions.
#'
#'   @note The dataset must contain the same headers (in the same order) that the \code{data}
#'         used to find the system with fugeR.run.
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
#'      prediction <- fugeR.predict( fis, In )
#'   }
#'
#' @seealso  \code{\link{fugeR.run}}
#' 
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.predict <- function( fuzzySystem=NULL, 
                           dataset = NULL) {
            
    #Fuzzify and apply AND (MIN) operator.
    #Return a list for each rule of length equal to the number of case.
    fuzzifiedValue <- fugeR.fuzzify(fuzzySystem, dataset)
    
    #Defuzzify and make prediction.
    #Return a lsit for each out of length equal to then number of case.
    prediction <- fugeR.defuzzify(fuzzySystem, fuzzifiedValue)
    
    return(data.frame(prediction))
} 
