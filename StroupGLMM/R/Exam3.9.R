#' @title Example 3.9 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-118)
#' @name   Exam3.9
#' @docType data
#' @keywords datasets
#' @description Exam3.9 used to differentiate conditional and marginal binomial models with and without interaction for S2 variable.
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
#' @importFrom MASS glmmPQL
#' @importFrom nlme lme
#' 
#' @examples
#' #-----------------------------------------------------------------------------------
#' ## Binomial conditional GLMM without interaction, logit link
#' #-----------------------------------------------------------------------------------
#' library(MASS)
#' DataSet3.2$trt <- factor( x  =  DataSet3.2$trt )
#' DataSet3.2$loc <- factor( x  =  DataSet3.2$loc )
#' Exam3.9.fm1   <-
#'   glmmPQL(
#'       fixed    =  S2/Nbin~trt
#'     , random   = ~1|loc
#'     , family   =  quasibinomial(link = "logit")
#'     , data     =  DataSet3.2
#'   # , weights
#'   # , control
#'     , niter    = 10
#'     , verbose  = TRUE
#'     # , ...
#'   )
#' summary(Exam3.9.fm1)
#' 
#' #-------------------------------------------------------------
#' ##  treatment means
#' #-------------------------------------------------------------
#' library(lsmeans)
#' (Lsm3.9fm1   <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.9.fm1
#'     , specs   = "trt"
#'     , link=TRUE
#'     # , ...
#'   )
#' )

#' ##--- Normal Approximation
#' library(nlme)
#' Exam3.9fm2 <-
#'   lme(
#'       fixed       = S2/Nbin~trt
#'     , data        = DataSet3.2
#'     , random      = ~1|loc
#'     , weights     = NULL
#'   # , subset
#'     , method      = "REML" #c("REML", "ML")
#'     , na.action   = na.fail
#'   # , control     = list()
#'     , contrasts   = NULL
#'     , keep.data   = TRUE
#'   )
#' (Lsm3.9fm2    <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.9fm2
#'     , specs   = "trt"
#'     # , ...
#'   )
#' )
#' 
#' ##---Binomial GLMM with interaction
#' Exam3.9fm3   <-
#'   glmmPQL(
#'       fixed       =  S2/Nbin~trt
#'     , random      = ~1|trt/loc
#'     , family      =  quasibinomial(link = "logit")
#'     , data        =  DataSet3.2
#'   # , weights
#'   # , control
#'     , niter = 10
#'     , verbose = TRUE
#'     # , ...
#'   )
#' summary(Exam3.9fm3)
#' (Lsm3.9fm3    <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.9fm3
#'     , specs   = "trt"
#'   # , ...
#'   )
#' )
#' 
#' ##---Binomial Marginal GLMM(assuming compound symmetry)
#' Exam3.9fm4   <-
#'   glmmPQL(
#'       fixed       =  S2/Nbin~trt
#'     , random      = ~1|loc
#'     , family      =  quasibinomial(link = "logit")
#'     , data        =  DataSet3.2
#'     , correlation =  corCompSymm(form=~1|loc)
#'   # , weights
#'   # , control
#'     , niter       = 10
#'     , verbose     = TRUE
#'   # , ...
#'   )
#' summary(Exam3.9fm4)
#' (Lsm3.9fm4  <-
#'   lsmeans::lsmeans(
#'       object  = Exam3.9fm4
#'     , specs   = "trt"
#'     # , ...
#'   )
#' )
NULL
