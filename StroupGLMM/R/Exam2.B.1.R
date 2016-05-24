#' @title Example 2.B.1 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-53)
#' @name   Exam2.B.1
#' @docType data
#' @keywords datasets
#' @description Exam2.B.1 is used to visualize the effect of lm model statement with Gaussian data and their design matrix
#' @author \enumerate{
#'          \item  Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Adeela Munawar (\email{adeela.uaf@@gmail.com})
#'          }
#' @references \enumerate{
#' \item Stroup, W. W. (2012).
#'      \emph{Generalized Linear Mixed Models: Modern Concepts, Methods and Applications}.
#'        CRC Press.
#'  }
#' @seealso
#'    \code{\link{Table1.1}}
#'    
#' @importFrom stats lm summary.lm
#' 
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## Linear Model  discussed in Example 2.B.1 using simple regression data of Table1.1
#' #-----------------------------------------------------------------------------------
#' data(Table1.1)
#' Exam2.B.1.lm1 <-
#'   lm(
#'       formula     = y~x
#'     , data        = Table1.1
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'     , method      = "qr"
#'     , model       = TRUE
#'     #  , x           = FALSE
#'     #  , y           = FALSE
#'     , qr          = TRUE
#'     , singular.ok = TRUE
#'     , contrasts   = NULL
#'     #  , offset
#'     #  , ...
#'   )
#' summary(Exam2.B.1.lm1)
#' DesignMatrix.lm1 <-
#'   model.matrix (
#'     object = Exam2.B.1.lm1
#'   )
#' DesignMatrix.lm1
NULL