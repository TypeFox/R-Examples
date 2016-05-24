#' Adaptive Neuro Fuzzy Inference System in R
#'
#' The package implements ANFIS Type 3 Takagi and Sugeno's fuzzy if-then rule 
#' network. This package includes the new following features: 
#' \enumerate{
#'   \item Membership Functions (MF) flexible framework:
#'    \itemize{
#'      \item  Flexible user-defined membership functions(MF) extensible class.
#'      \item  Independent number of (MF) for each input.   
#'      \item  Different MF types, if required, for each input. 
#'    }
#'   \item Type 3 Takagi and Sugeno's fuzzy if-then rule 
#'   \item Full Rule combinations, e.g. 2 inputs 2 membership functions this 
#'     means that 4 fuzzy rules will be created.
#'   \item Different learning strategies:
#'    \describe{
#'      \item{trainHybridJangOffLine}{Hybrid learning, i.e. Descent Gradient 
#'        for precedents and Least Squares Estimation for consequents.}
#'      \item{trainHybridJangOnLine}{on-line version with hybrid learning.}
#'      \item{trainHybridOffLine}{Adaptive learning coefficient and momentum 
#'        term.}
#'     }  
#'   \item Multiple outputs support, i.e., the same input partition can be used
#'    to predict more than one output variable.
#' }
#'
#' @docType package
#' @name Anfis-package
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
#' @keywords ANFIS membership fuzzy
#' @references 
#' \enumerate{
#'  \item Jang, J. S. (1993). ANFIS: adaptive-network-based fuzzy inference 
#'  system. Systems, Man and Cybernetics, IEEE Transactions on, 23(3), 665-685.
#' }
NULL