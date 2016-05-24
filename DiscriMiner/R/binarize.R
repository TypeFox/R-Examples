#' Binarize a data frame into a super-indicator matrix
#' 
#' Convert a data frame with factors into a super-indicator matrix (a.k.a.
#' complete disjunctive table from the french \emph{tableau disjonctive
#' complete})
#' 
#' 
#' @param variables data frame with categorical variables (coded as factors)
#' @return A super-indicator matrix of binary data
#' @author Gaston Sanchez
#' @seealso \code{\link{easyMCA}}
#' @references Lebart L., Piron M., Morineau A. (2006) \emph{Statistique
#' Exploratoire Multidimensionnelle}. Dunod, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load insurance cars dataset
#'   data(insurance)
#' 
#'   # super-indicator matrix of binary data
#'   bin_insure = binarize(insurance[,-1])
#'   head(bin_insure)
#'   }
#' 
binarize <- 
function(variables)
{
  # binary super-indicator matrix (aka Complete Disjunctive Table)
  # variables: matrix or data.frame with explanatory variables
  
  # make sure variables is a data frame with factors
  fac_check = sapply(variables, class)
  if (!is.data.frame(variables) && any(fac_check != "factor"))
    stop("\n'variables' must be a data frame with factors")
  # no missing values allowed
  if (length(complete.cases(variables)) != nrow(variables))
    stop("\nSorry, no missing values allowed in 'variables'")    
  
  # build super-indicator matrix Z
  Z = my_tdc(variables)  
  Z
}
