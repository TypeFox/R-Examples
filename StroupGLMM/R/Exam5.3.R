#' @title Example 5.3 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-172)
#' @name   Exam5.3
#' @docType data
#' @keywords datasets
#' @description Exam5.3 Inference using empirical standard error with different Bias connection
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
#' @importFrom  lme4 lmer
#' 
#' 
#' @examples
#' 
#' data(DataSet4.1)
#' DataSet4.1$trt   <- factor(x = DataSet4.1$trt)
#' DataSet4.1$block <- factor( x = DataSet4.1$block)
#' 
#' ##---REML estimates on page 172
#' library(lme4)
#' # library(lmerTest)
#' Exam5.3REML      <-
#'   lmer(
#'          formula    = y ~ trt + (1|block)
#'        , data       = DataSet4.1
#'        , REML       = TRUE
#'     #  , control    = lmerControl()
#'        , start      = NULL
#'     #  , verbose    = 0L
#'     #  , subset
#'     #  , weights
#'     #  , na.action
#'     #  , offset
#'        , contrasts  = NULL
#'        , devFunOnly = FALSE
#'     #  , ...
#'   )
#' ##---Standard Error Type "Model Based" with no Bias Connection
#' AnovaExam5.3REML  <- anova( object = Exam5.3REML )
#' AnovaExam5.3REML
#' 
#' ##---Standard Error Type "Model Based" with "Kenward-Roger approximation" Bias Connection
#' # library(pbkrtest)
#' anova( object = Exam5.3REML, ddf  = "Kenward-Roger")
#' 
#' ##---ML estimates on page 172
#' Exam5.3ML      <-
#'   lmer(
#'         formula      = y ~ trt + ( 1|block )
#'        , data        = DataSet4.1
#'        , REML        = FALSE
#'    #   , control     = lmerControl()
#'        , start       = NULL
#'    #   , verbose     = 0L
#'    #   , subset
#'    #   , weights
#'    #   , na.action
#'    #   , offset
#'        , contrasts   = NULL
#'        , devFunOnly  = FALSE
#'    #   , ...
#'   )
#'   
#' ##---Standard Error Type "Model Based" with no Bias Connection
#' AnovaExam5.3ML  <- anova( object = Exam5.3ML )
#' AnovaExam5.3ML
#' 
#' ##---Standard Error Type "Model Based" with "Kenward-Roger approximation" Bias Connection
#' anova( object = Exam5.3ML, ddf = "Kenward-Roger")
NULL