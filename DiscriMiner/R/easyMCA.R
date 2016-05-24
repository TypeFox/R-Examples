#' Multiple Correspondence Analysis
#' 
#' Performs a basic Multiple Correspondence Analysis (MCA)
#' 
#' 
#' @param variables data frame with categorical variables (coded as factors)
#' @return An object of class \code{"qualmca"}, basically a list with the
#' following elements:
#' @return \item{values}{table with eigenvalues}
#' @return \item{coefficients}{coefficients of factorial axes}
#' @return \item{components}{factor coordinates}
#' @author Gaston Sanchez
#' @seealso \code{\link{disqual}}, \code{\link{binarize}}
#' @references Lebart L., Piron M., Morineau A. (2006) \emph{Statistique
#' Exploratoire Multidimensionnelle}. Dunod, Paris.
#' 
#' Saporta G. (2006) \emph{Probabilites, analyse des donnees et statistique}.
#' Editions Technip, Paris.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load insurance wines dataset
#'   data(insurance)
#' 
#'   # multiple correspondence analysis
#'   mca1 = easyMCA(insurance[,-1])
#'   mca1
#'   }
#' 
easyMCA <- 
function(variables)
{
  # Perform multiple correspondence analysis
  # X: data frame with categorical variables as factors
  
  # check input
  fac_check = sapply(variables, class)
  if (!is.data.frame(variables) && any(fac_check != "factor"))
    stop("\nA data frame with factors was expected")
  # check for missing values
  if (length(complete.cases(variables)) != nrow(variables))
    stop("\nSorry: no missing values allowed")
  
  # apply MCA
  res = my_mca(variables)
  res
}
