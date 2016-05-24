#' @title Example 3.5 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-85)
#' @name   Exam3.5
#' @docType data
#' @keywords datasets
#' @description Exam3.5 fixed location, factorial treatment structure, Gaussian response
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
#'    \code{\link{DataSet3.2}}
#'    
#' @importFrom lsmeans lsmeans contrast
#' @importFrom phia testFactors
#' @importFrom stats glm summary.glm
#' @importFrom mutoss sidak
#' 
#' @examples
#' data(DataSet3.2)
#' DataSet3.2$A <- factor(x = DataSet3.2$A)
#' DataSet3.2$B <- factor(x = DataSet3.2$B)
#' DataSet3.2$loc <- factor(x = DataSet3.2$loc, level = c(8, 1, 2, 3, 4, 5, 6, 7))
#' Exam3.5.lm <-
#'   lm(
#'          formula     =  Y~ A + B +loc
#'        , data        =  DataSet3.2
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'        , method      =  "qr"
#'        , model       =  TRUE
#'     #  , x           =  FALSE
#'     #  , y           =  FALSE
#'        , qr          =  TRUE
#'        , singular.ok =  TRUE
#'        , contrasts   =  NULL
#'     #  , offset
#'     #  , ...
#'   )
#'   
#' ##---a0 marginal mean
#' list3.5.a0 <- list(B = c("0" = 1,"1" = 0) )
#' library(phia)
#' Test3.5.a0 <-
#'   summary(testFactors(
#'       model  =  Exam3.5.lm
#'     , levels =  list3.5.a0)
#'     )
#'     
#' ##---b0 marginal mean
#' list3.5.b0 <- list(B = c("0" = 1,"1" = 0) )
#' Test3.5.b0 <-
#' summary(testFactors(
#'     model  =  Exam3.5.lm
#'   , levels =  list3.5.b0)
#'   )
#'   
#' ##---Simple effect of A on B0
#' Test3.5.AB0 <-
#'   summary(testInteractions(
#'       model  =  Exam3.5.lm
#'     , custom =  list3.5.b0
#'     , across =  "B")
#'     )
#'     
#' ##---Simple effect of B on A0
#' Test3.5.BA0 <-
#'   summary(testInteractions(
#'       model  =  Exam3.5.lm
#'     , custom =  list3.5.a0
#'     , across =  "A")
#'     )
#'     
#' ##---Simple Effect of A over B
#' (SimpleEffect3.5.AB <-
#'   summary(testInteractions(
#'       model  =  Exam3.5.lm
#'     , fixed  =   "A"
#'     , across =  "B")
#'     )
#'     )
#'     
#' ##---Simple Effect of B over A
#' (SimpleEffect3.5.BA <-
#'   summary(testInteractions(
#'       model  =  Exam3.5.lm
#'     , fixed  =   "B"
#'     , across =  "A")
#'     )
#' )
#' #-------------------------------------------------------------
#' ## Individula least squares treatment means
#' #-------------------------------------------------------------
#' (Lsm3.5 <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.5.lm
#'     , specs   = ~A*B
#'     # , ...
#'   )
#' )
NULL
