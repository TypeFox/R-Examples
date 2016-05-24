#' @title Example 2.B.5 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-57)
#' @name   Exam2.B.5
#' @docType data
#' @keywords datasets

#' @description Exam2.B.5 is related to multi batch regression data assuming different forms of linear models.
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
#' @importFrom stats lm summary.lm model.matrix  
#'    
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## Nested Model with no intercept
#' #-----------------------------------------------------------------------------------
#' data(Table1.2)
#' Table1.2$Batch <- factor(x = Table1.2$Batch)
#' Exam2.B.5.lm1 <-
#'   lm(
#'          formula     = Y~0+Batch+ Batch/X
#'        , data        = Table1.2
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'        , method      = "qr"
#'        , model       = TRUE
#'     #  , x           = FALSE
#'     #  , y           = FALSE
#'        , qr          = TRUE
#'        , singular.ok = TRUE
#'        , contrasts   = NULL
#'     #  , offset
#'     #  , ...
#'   )
#' DesignMatrix.lm1 <- model.matrix (object = Exam2.B.5.lm1)
#' DesignMatrix.lm1
#' #-----------------------------------------------------------------------------------
#' ## Interaction Model with intercept
#' #-----------------------------------------------------------------------------------
#' Exam2.B.5.lm2 <-
#'   lm(
#'          formula     = Y~Batch +X+ Batch*X
#'        , data        = Table1.2
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'        , method      = "qr"
#'        , model       = TRUE
#'     #  , x           = FALSE
#'     #  , y           = FALSE
#'        , qr          = TRUE
#'        , singular.ok = TRUE
#'        , contrasts   = NULL
#'     #  , offset
#'     #  , ...
#'   )
#' DesignMatrix.lm2 <-   model.matrix (object = Exam2.B.5.lm2)
#' DesignMatrix.lm2
#' #-----------------------------------------------------------------------------------
#' ## Interaction Model with no intercept
#' #-----------------------------------------------------------------------------------
#' Exam2.B.5.lm3 <-
#'   lm(
#'          formula     = Y~0 + Batch + Batch*X
#'        , data        = Table1.2
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'        , method      = "qr"
#'        , model       = TRUE
#'     #  , x           = FALSE
#'     #  , y           = FALSE
#'        , qr          = TRUE
#'        , singular.ok = TRUE
#'        , contrasts   = NULL
#'     #  , offset
#'     #  , ...
#'   )
#' DesignMatrix.lm3 <-   model.matrix(object = Exam2.B.5.lm3)
#' #-----------------------------------------------------------------------------------
#' ## Interaction Model with intercept  but omitting X term as main effect
#' #-----------------------------------------------------------------------------------
#' Exam2.B.5.lm4 <-
#'   lm(
#'          formula     = Y~Batch + Batch*X
#'        , data        = Table1.2
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'        , method      = "qr"
#'        , model       = TRUE
#'     #  , x           = FALSE
#'     #  , y           = FALSE
#'        , qr          = TRUE
#'        , singular.ok = TRUE
#'        , contrasts   = NULL
#'     #  , offset
#'     #  , ...
#'   )
#' DesignMatrix.lm4 <-   model.matrix(object = Exam2.B.5.lm4)
#' DesignMatrix.lm4
NULL