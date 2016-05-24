#' @title Example 2.B.4 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-56)
#' @name   Exam2.B.4
#' @docType data
#' @keywords datasets
#' @description Exam2.B.4 is used to illustrate one way treatment design with Binomial observations.
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
#'    \code{\link{DataExam2.B.4}}
#'    
#' @importFrom stats glm summary.glm
#'        
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## logit Model  discussed in Example 2.B.2 using DataExam2.B.4
#' ## Default link is logit
#' ## using fmaily=binomial gives warning message of no-integer successes
#' #-----------------------------------------------------------------------------------
#' data(DataExam2.B.4)
#' DataExam2.B.4$trt <- factor(x =  DataExam2.B.4$trt)
#' Exam2.B.4glm <-
#'   glm(
#'          formula    =  Yij/Nij~trt
#'        , family     =  quasibinomial(link = "probit")
#'        , data       =  DataExam2.B.4
#'        , weights    =  NULL
#'     #  , subset
#'     #  , na.action
#'        , start      =  NULL
#'     #  , etastart
#'     #  , mustart
#'     #  , offset
#'     #  , control    =  list(...)
#'     #  , model      =  TRUE
#'        , method     =  "glm.fit"
#'     #  , x          =  FALSE
#'     #  , y          =  TRUE
#'        , contrasts  =  NULL
#'     #  , ...
#'   )
#' summary(Exam2.B.4glm)
NULL