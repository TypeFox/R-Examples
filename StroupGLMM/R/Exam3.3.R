#' @title Example 3.3 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-77)
#' @name   Exam3.3
#' @docType data
#' @keywords datasets
#' @description Exam3.3 use RCBD data with fixed location effect and different forms of estimable functions are shown in this example.
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
#' #-----------------------------------------------------------------------------------
#' ## linear model for Gaussian data
#' #-----------------------------------------------------------------------------------
#' data(DataSet3.2)
#' DataSet3.2$trt <- factor(x = DataSet3.2$trt, level = c(3,0,1,2))
#' DataSet3.2$loc <- factor(x = DataSet3.2$loc, level = c(8, 1, 2, 3, 4, 5, 6, 7))
#' Exam3.3.lm1 <-
#'   lm(
#'          formula     = Y~ trt+loc
#'        , data        = DataSet3.2
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
#' summary( Exam3.3.lm1 )
#' #-------------------------------------------------------------
#' ## Individula least squares treatment means
#' #-------------------------------------------------------------
#' library(lsmeans)
#' (Lsm3.3    <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.3.lm1
#'     , specs   = "trt"
#'     # , ...
#'   )
#' )
#' #---------------------------------------------------
#' ## Pairwise treatment means estimate
#' #---------------------------------------------------
#' contrast( object = Lsm3.3 , method = "pairwise")
#' #---------------------------------------------------
#' ## Repairwise treatment means estimate
#' #---------------------------------------------------
#' ## contrast( object = Lsm3.3 , method = "repairwise")
#' #-------------------------------------------------------
#' ## LSM Trt0 (This term is used in Walter Stroups' book)
#' #-------------------------------------------------------
#' library(phia)
#' list3.3.1 <- list(trt=c("0" = 1 ) )
#' Test3.3.1 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels =  list3.3.1)
#'   )
#' #-------------------------------------------------------
#' ## LSM Trt0 alt(This term is used in Walter Stroups' book)
#' #-------------------------------------------------------
#' list3.3.2 <-
#'   list(trt=c("0" = 1 )
#'        , loc=c("1" = 0,"2" = 0,"3" = 0,"4" = 0,"5" = 0,"6" = 0,"7" = 0)
#'   )
#' Test3.3.2 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels =  list3.3.2)
#'   )
#' #-------------------------------------------------------
#' ##  Trt0 Vs Trt1
#' #-------------------------------------------------------
#' list3.3.3 <- list(trt=c("0" = 1,"1" = -1))
#' Test3.3.3 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels =  list3.3.3)
#'   )
#' #-------------------------------------------------------
#' ##  average Trt0+1
#' #-------------------------------------------------------
#' list3.3.4 <- list(trt=c("0" = 0.5 , "1" = 0.5))
#' Test3.3.4 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels =  list3.3.4)
#'   )
#' #-------------------------------------------------------
#' ##  average Trt0+2+3
#' #-------------------------------------------------------
#' list3.3.5 <- list(trt=c("0" = 0.33333,"2" = 0.33333,"3" = 0.33333))
#' Test3.3.5 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.5)
#'   )
#' #-------------------------------------------------------
#' ##  Trt 2 Vs 3 difference
#' #-------------------------------------------------------
#' list3.3.6 <- list(trt=c("2" = 1,"3" = -1))
#' Test3.3.6 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.6)
#'   )
#' #-------------------------------------------------------
#' ##  Trt 1 Vs 2 difference
#' #-------------------------------------------------------
#' list3.3.7 <- list(trt=c("1" = 1,"2" = -1))
#' Test3.3.7 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.7)
#'   )
#' #-------------------------------------------------------
#' ##  Trt 1 Vs 3 difference
#' #-------------------------------------------------------
#' list3.3.8 <- list(trt=c("1" = 1,"3" = -1))
#' Test3.3.8 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.8)
#'   )
#' #-------------------------------------------------------
#' ##  Average trt0+1  vs Average Trt2+3
#' #-------------------------------------------------------
#' list3.3.9 <-  list(trt=c("0" = 0.5,"1" = 0.5,"2" = -0.5,"3" = -0.5))
#' Test3.3.9 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.9)
#'   )
#' #-------------------------------------------------------
#' ##  Trt1  vs Average Trt0+1+2
#' #-------------------------------------------------------
#' list3.3.10 <- list(trt=c("0" = 0.33333,"1" = -1,"2" = 0.33333,"3" = 0.33333))
#' Test3.3.10 <-
#' summary(testFactors(
#'     model  =  Exam3.3.lm1
#'   , levels = list3.3.10)
#'   )
#' #-------------------------------------------------------
#' ## Sidak Multiplicity adjustment for p-values
#' #-------------------------------------------------------
#' library(mutoss)
#' PValues3.3 <-
#'   c(
#'     Test3.3.3[[7]][1, 4]
#'   , Test3.3.6[[7]][1, 4]
#'   , Test3.3.7[[7]][1, 4]
#'   , Test3.3.8[[7]][1, 4]
#'   , Test3.3.9[[7]][1, 4]
#'   , Test3.3.10[[7]][1, 4]
#'    )
#'  AdjPValues3.3 <- sidak(PValues3.3)
NULL