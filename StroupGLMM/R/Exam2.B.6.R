#' @title Example 2.B.6 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-58)
#' @name   Exam2.B.6
#' @docType data
#' @keywords datasets
#' @description Exam2.B.6 is related to multi batch regression data assuming different forms of linear models keeping batch effect random.
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
#'    \code{\link{Table1.2}}
#' 
#' @importFrom nlme lme
#' 
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## Nested Model with no intercept
#' #-----------------------------------------------------------------------------------
#' data(Table1.2)
#' library(nlme)
#' Table1.2$Batch <- factor(x = Table1.2$Batch)
#' Exam2.B.6fm1 <-
#'   lme(
#'       fixed       = Y~X
#'     , data        = Table1.2
#'     , random      = list(Batch = pdDiag(~1), X = pdDiag(~1))
#'     , correlation = NULL
#'     , weights     = NULL
#'   # , subset
#'     , method      = "REML" #c("REML", "ML")
#'     , na.action   = na.fail
#'   # , control     = list()
#'     , contrasts   = NULL
#'     , keep.data   = TRUE
#'   )
NULL
