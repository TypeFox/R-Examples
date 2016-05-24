#' @title Example 2.B.7 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-60)
#' @name   Exam2.B.7
#' @docType data
#' @keywords datasets
#' @description Exam2.B.7 is related to multi batch regression data assuming different forms of linear models with factorial experiment.
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
#'    \code{\link{DataExam2.B.7}}
#'    
#' @importFrom stats lm summary.lm model.matrix lm.fit coef
#'    
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## Classical main effects and Interaction Model
#' #-----------------------------------------------------------------------------------
#' data(DataExam2.B.7)
#' DataExam2.B.7$a <- factor(x = DataExam2.B.7$a)
#' DataExam2.B.7$b <- factor(x = DataExam2.B.7$b)
#' Exam2.B.7.lm1 <-
#'   lm(
#'         formula     = y~ a + b + a*b
#'       , data        = DataExam2.B.7
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
#' #-----------------------------------------------------------------------------------
#' ## One way treatment effects model
#' #-----------------------------------------------------------------------------------
#' DesignMatrix.lm1 <- model.matrix (object = Exam2.B.7.lm1)
#' DesignMatrix2.B.7.2 <- DesignMatrix.lm1[,!colnames(DesignMatrix.lm1) %in% c("a2","b")]
#' lmfit2 <-
#'   lm.fit(
#'       x           = DesignMatrix2.B.7.2
#'     , y           = DataExam2.B.7$y
#'     , offset      = NULL
#'     , method      = "qr"
#'     , tol         = 1e-07
#'     , singular.ok = TRUE
#'     # , ...
#'   )
#' Coefficientslmfit2 <- coef( object = lmfit2)
#' #-----------------------------------------------------------------------------------
#' ## One way treatment effects model without intercept
#' #-----------------------------------------------------------------------------------
#' DesignMatrix2.B.7.3    <-
#'   as.matrix(DesignMatrix.lm1[,!colnames(DesignMatrix.lm1) %in% c("(Intercept)","a2","b")])
#'   
#' lmfit3 <-
#'   lm.fit(
#'       x           = DesignMatrix2.B.7.3
#'     , y           = DataExam2.B.7$y
#'     , offset      = NULL
#'     , method      = "qr"
#'     , tol         = 1e-07
#'     , singular.ok = TRUE
#'     # , ...
#'   )
#' Coefficientslmfit3 <- coef( object = lmfit3)
#' 
#' #-----------------------------------------------------------------------------------
#' ## Nested Model (both models give the same result)
#' #-----------------------------------------------------------------------------------
#' Exam2.B.7.lm4 <-
#'   lm(
#'          formula     = y~ a + a/b
#'        , data        = DataExam2.B.7
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
#' summary(Exam2.B.7.lm4)
#' 
#' Exam2.B.7.lm4 <-
#'   lm(
#'          formula     = y~ a + a*b
#'        , data        = DataExam2.B.7
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
#' summary(Exam2.B.7.lm4)
NULL
