#' @title Example 2.B.2 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-54)
#' @name   Exam2.B.2
#' @docType data
#' @keywords datasets
#' @description Exam2.B.2 is used to visualize the effect of glm model statement with binomial data with logit and probit links.
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
#'    \code{\link{DataExam2.B.2}}
#'    
#' @importFrom stats glm summary.glm
#' 
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## probitit Model  discussed in Example 2.B.2 using DataExam2.B.2
#' ## Default link is logit
#' ## using fmaily=binomial gives warning message of no-integer successes
#' #-----------------------------------------------------------------------------------
#' data(DataExam2.B.2)
#' Exam2.B.2glm <-
#'   glm(
#'       formula    =  y/n~x
#'     , family     =  quasibinomial(link = "probit")
#'     , data       =  DataExam2.B.2
#'     , weights    =  NULL
#'  #  , subset
#'  #  , na.action
#'     , start      =  NULL
#'  #  , etastart
#'  #  , mustart
#'  #  , offset
#'  #  , control    =  list(...)
#'  #  , model      =  TRUE
#'     , method     =  "glm.fit"
#'  #  , x          =  FALSE
#'  #  , y          =  TRUE
#'     , contrasts  =  NULL
#'     #  , ...
#'   )
#' summary(Exam2.B.2glm)
NULL
