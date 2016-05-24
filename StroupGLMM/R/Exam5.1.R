#' @title Example 5.1 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-163)
#' @name   Exam5.1
#' @docType data
#' @keywords datasets
#' @description Exam5.1 is used to show polynomial multiple regression with binomial response
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
#'    \code{\link{DataSet5.1}}
#'    
#' @importFrom aod wald.test
#' @importFrom stats lm summary.lm glm summary.glm anova
#' 
#' @examples
#' 
#' ##---Sequential Fit of the logit Model
#' Exam5.1.glm.1 <-
#'   glm(
#'       formula    =  F/N~ X
#'     , family     =  quasibinomial(link = "logit")
#'     , data       =  DataSet5.1
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
#'  #  , ...
#'   )
#' summary(Exam5.1.glm.1)
#' 
#' ## confint.default()   produce Wald Confidence interval as SAS produces
#' ##---Likelihood Ratio test for Model 1
#' (LRExam5.1.glm.1  <-
#'   anova(
#'       object =  Exam5.1.glm.1
#'     , test   = "Chisq")
#' )
#' 
#' library(aod)
#' WaldExam5.1.glm.1 <-
#'   wald.test(
#'       Sigma   = vcov(object=Exam5.1.glm.1)
#'     , b       = coef(object=Exam5.1.glm.1)
#'     , Terms   = 2
#'     , L       = NULL
#'     , H0      = NULL
#'     , df      = NULL
#'     , verbose = FALSE
#'   )
#' 
#' ##---Sequential Fit of the logit Model quadratic terms involved
#' Exam5.1.glm.2 <-
#'   glm(
#'       formula    =  F/N~ X + I(X^2)
#'     , family     =  quasibinomial(link = "logit")
#'     , data       =  DataSet5.1
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
#'  #  , ...
#'   )
#' summary( Exam5.1.glm.2 )
#' 
#' ##---Likelihood Ratio test for Model Exam5.1.glm.2
#' (LRExam5.1.glm.2   <-
#'   anova(
#'       object =  Exam5.1.glm.2
#'     , test   = "Chisq")
#' )
#' 
#' WaldExam5.1.glm.2 <-
#'   wald.test(
#'       Sigma   = vcov(object=Exam5.1.glm.2)
#'     , b       = coef(object=Exam5.1.glm.2)
#'     , Terms   = 3
#'     , L       = NULL
#'     , H0      = NULL
#'     , df      = NULL
#'     , verbose = FALSE
#'   )
#'   
#' ##---Sequential Fit of the logit Model 5th power terms involved
#' Exam5.1.glm.3 <-
#'   glm(
#'       formula    =  F/N~ X + I(X^2) + I(X^3) + I(X^4) + I(X^5)
#'     , family     =  quasibinomial(link = "logit")
#'     , data       =  DataSet5.1
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
#'  #  , ...
#'   )
#' summary(Exam5.1.glm.3)
#' 
#' ## confint.default()   produce Wald Confidence interval as SAS produces
#' ##---Likelihood Ratio test for Model 1
#' (LRExam5.1.glm.3   <-
#'   anova(
#'       object =  Exam5.1.glm.3
#'     , test   = "Chisq")
#' )
#' 
#' WaldExam5.1.glm.3 <-
#'   wald.test(
#'       Sigma   = vcov(object=Exam5.1.glm.3)
#'     , b       = coef(object=Exam5.1.glm.3)
#'     , Terms   = 6
#'     , L       = NULL
#'     , H0      = NULL
#'     , df      = NULL
#'     , verbose = FALSE
#'   )
NULL
