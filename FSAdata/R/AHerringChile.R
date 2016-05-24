#' @title Ages and lengths of Araucanian Herring from Chilean waters.
#' 
#' @description Ages and lengths of Araucanian Herring (\emph{Strangomera bentincki}) from Chilean waters.
#' 
#' @name AHerringChile
#' 
#' @docType data
#' 
#' @format A data frame with the following 2 variables:
#'  \describe{
#'    \item{age}{Age in years.}
#'    \item{len}{Total length (to nearest 0.5 cm).}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item Seasonal Growth
#'    \item von Bertalanffy 
#'    \item Somers model
#'  }
#'  
#' @concept 'Seasonal Growth' 'von Bertalanffy' Somers
#' 
#' @source From figure 9 of Cubillos, L.A., D.F. Arcosa, D.A. Bucareya, M.T. Canalesa.  2001.  Seasonal growth of small pelagic fish off Talcahuano, Chile (37S, 73W): a consequence of their reproductive strategy to seasonal upwelling?  Aquatic Living Resources, 14:115-124.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(AHerringChile)
#' str(AHerringChile)
#' head(AHerringChile)
#' plot(len~age,data=AHerringChile)
#' 
NULL
