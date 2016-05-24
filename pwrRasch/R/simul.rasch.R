#' Simulate data according to the Rasch model
#' 
#' This function simulates data according to the Rasch model
#' based on user-specified item and person parameters.
#' 
#' If persons is an integer value, the corresponding parameter vector 
#' is drawn from N(0, 1.5). If items is an integer value, the corresponding parameter vector
#' is equally spaced between [-3, 3]. Note that item parameters need to be normalized to sum-0. 
#' This precondition can be overruled using argument \code{sum0 = FALSE}.
#'   
#' @param persons     Either a vector of specified person parameters or an integer indicating the number of persons.
#' @param items       Either a vector of specified item parameters or an integer indicating the number of items. 
#' @param sum0        If \code{TRUE}, specified item parameters need to be normalized to sum-0.
#'  
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#'  
#' @seealso 
#' \code{\link{aov.rasch}}, \code{\link{pwr.rasch}}
#'
#' @references
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model 
#' calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model.
#' \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' @return 
#' Returns a 0-1 matrix according to the Rasch model.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # simulate Rasch model based data
#' # 100 persons, 20 items,
#' # person parameter drawn from a normal distribution: N(0,1.5)
#' # item parameters equally spaced between [-3, 3]
#' simul.rasch(100, items = 20)
#'
#' # simulate Rasch model based data
#' # 100 persons, 17 items
#' # person parameter drawn from a uniform distribution: U[-4, 4]
#' # item parameters: [-4.0, -3.5, -3.0, ... , 3.0, 3.5, 4.0]
#' simul.rasch(runif(100, -4, 4), items = seq(-4, 4, by = 0.5))
#' }
simul.rasch <- function(persons, items, sum0 = TRUE) {

  #------------------------------------------------------------------------------------------------------#
  # Input Check
  
  # sum of item parameters = 0
  if (length(items) == 1) {
    
    if (sum0 == TRUE & round(sum(seq(-3, 3, length.out = items)), digits = 3) != 0) {
      
      stop("Item pararameters are not normalized to sum-0")
      
    } 
  
  } else {
    
    if (sum0 == TRUE & sum(round(items, digits = 3)) != 0) {
      
      stop("Item pararameters are not normalized to sum-0")
      
    } 
    
  } 
  
  #------------------------------------------------------------------------------------------------------#
  
  # item parameters
  if (length(items) == 1) {
    
    diff <- seq(-3, 3, length.out = items)
    n.items <- items
    
  } else {
    
    diff <- items
    n.items <- length(items)
    
  }
  
  # person parameters
  if (length(persons) == 1) {
    
    ability <- rnorm(persons, sd = 1.5)
    n.persons <- persons
    
  } else {
    
    ability <- persons
    n.persons <- length(persons)
    
  }
  
  fsmat <- outer(ability, diff, "-")
  psolve <- exp(fsmat) / (1 + exp(fsmat))
  
  resmat <- (matrix(runif(n.items * n.persons), n.persons, n.items) < psolve) * 1
  
  return(resmat)
  
}