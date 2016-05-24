#' @title Example 4.1 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-138)
#' @name   Exam4.1
#' @docType data
#' @keywords datasets
#' @description Exam4.1 REML vs ML criterion is used keeping block effects random
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
#'    \code{\link{DataSet4.1}}
#'    
#' @importFrom lme4 lmer
#' 
#' @examples
#' 
#' DataSet4.1$trt   <- factor(x =  DataSet4.1$trt)
#' DataSet4.1$block <- factor(x =  DataSet4.1$block)
#'
#'##---REML estimates on page 138(article 4.4.3.3)
#'library(lme4)
#' Exam4.1REML  <-
#'   lmer(
#'       formula     = y~ trt +( 1|block )
#'     , data        = DataSet4.1
#'     , REML        = TRUE
#' #  , control     = lmerControl()
#'     , start       = NULL
#' #  , verbose     = 0L
#' #  , subset
#' #  , weights
#' #  , na.action
#' #  , offset
#'     , contrasts   = NULL
#'     , devFunOnly  = FALSE
#' #  , ...
#'   )
#'   
#' VarCompREML4.1  <-
#'   VarCorr(x     =   Exam4.1REML
#'           # , sigma = 1
#'           # , ...
#'   )
#' print(VarCompREML4.1, comp=c("Variance"))
#' 
#' ##---ML estimates on page 138(article 4.4.3.3)
#' Exam4.1ML  <-
#'   lmer(
#'         formula     = y ~ trt + (1|block)
#'        , data       = DataSet4.1
#'        , REML       = FALSE
#'    #  , control     = lmerControl()
#'        , start      = NULL
#'    #  , verbose     = 0L
#'    #  , subset
#'    #  , weights
#'    #  , na.action
#'    #  , offset
#'        , contrasts   = NULL
#'        , devFunOnly  = FALSE
#'    #  , ...
#'   )
#' VarCompML4.1  <-
#'   VarCorr(x     =    Exam4.1ML
#'           # , sigma = 1
#'           # , ...
#'   )
#' print(VarCompML4.1,comp=c("Variance"))
#' 
#' Exam4.1.lm <-
#'   lm(
#'       formula     = y~ trt + block
#'     , data        = DataSet4.1
#'  #  , subset
#'  #  , weights
#'  #  , na.action
#'     , method      = "qr"
#'     , model       = TRUE
#'  #  , x           = FALSE
#'  #  , y           = FALSE
#'     , qr          = TRUE
#'     , singular.ok = TRUE
#'     , contrasts   = NULL
#'  #  , offset
#'  #  , ...
#'   )
#' summary(anova(object = Exam4.1.lm))
NULL
