#' @title Example1.1 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-5)
#' @name   Exam1.1
#' @docType data
#' @keywords datasets
#' @description Exam1.1 is used for inspecting probability distribution and to define a plausible process through
#' linear models and generalized linear models.
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
#'    \code{\link{Table1.1}}
#'    
#' @importFrom ggplot2 ggplot
#' @importFrom stats lm summary.lm glm summary.glm cor
#' @importFrom survey svydesign svyglm
#' 
#' @examples
#' #-------------------------------------------------------------
#' ## Linear Model and results discussed in Article 1.2.1 after Table1.1
#' #-------------------------------------------------------------
#' data(Table1.1)
#' Exam1.1.lm1 <-
#'   lm(
#'       formula     =  y/Nx~x
#'     , data        =  Table1.1
#'     # , subset
#'     # , weights
#'     # , na.action
#'     , method      =  "qr"
#'     , model       =  TRUE
#'     , x           =  FALSE
#'     , y           =  FALSE
#'     , qr          =  TRUE
#'     , singular.ok =  TRUE
#'     , contrasts   =  NULL
#'     # , offset
#'     # , ...
#'   )
#' summary(Exam1.1.lm1 )
#' 
#' #-------------------------------------------------------------
#' ## GLM fitting with logit link (family=binomial)
#' #-------------------------------------------------------------
#' Exam1.1.glm1 <-
#'   glm(
#'       formula     =  y/Nx~x
#'     , family      =  binomial(link = "logit")
#'     , data        =  Table1.1
#'     , weights     =  NULL
#'     # , subset
#'     # , na.action
#'     , start       =  NULL
#'     # , etastart
#'     # , mustart
#'     # , offset
#'     # , control     =  list(...)
#'     # , model       =  TRUE
#'     , method      =  "glm.fit"
#'     , x           =  FALSE
#'     , y           =  TRUE
#'     , contrasts   =  NULL
#'     # , ...
#'   )
#' ## this glm() function gives warning message of non-integer success
#' summary(Exam1.1.glm1)
#' 
#' #-------------------------------------------------------------
#' ## GLM fitting with logit link (family=Quasibinomial)
#' #-------------------------------------------------------------
#' Exam1.1.glm2 <-
#'   glm(
#'       formula   =  y/Nx~x
#'     , family    =  quasibinomial(link = "logit")
#'     , data      =  Table1.1
#'     , weights   =  NULL
#'     # , subset
#'     # , na.action
#'     , start     =  NULL
#'     # , etastart
#'     # , mustart
#'     # , offset
#'     # , control   =  list(...)
#'     # , model     =  TRUE
#'     , method    =  "glm.fit"
#'     , x         =  FALSE
#'     , y         =  TRUE
#'     , contrasts =  NULL
#'     # , ...
#'   )
#' ## problem of "warning message of non-integer success" is overome by using quasibinomial family
#' summary(Exam1.1.glm2)
#' 
#' #-------------------------------------------------------------
#' ## GLM fitting with survey package(produces same result as using quasi binomial family in glm)
#' #-------------------------------------------------------------
#' library(survey)
#' design   <-
#'   svydesign(
#'       ids          =  ~1
#'     , probs        =  NULL
#'     , strata       =  NULL
#'     , variables    =  NULL
#'     , fpc          =  NULL
#'     , data         =  Table1.1
#'     # , nest         =  FALSE
#'     # , check.strata =  !nest
#'     , weights      =  NULL
#'     , pps          =  FALSE
#'     # , ...
#'   )
#' 
#' Exam1.1.svyglm  <-
#'   svyglm(
#'       formula  =  y/Nx~x
#'     , design   =  design
#'     # , ...
#'     , family   =  quasibinomial(link="logit")
#'   )
#' # summary(Exam1.1.svyglm)
#' 
#' #-------------------------------------------------------------
#' ## Figure 1.1
#' #-------------------------------------------------------------
#' Newdata     <-
#'   data.frame(
#'     Table1.1
#'     , LM       =  Exam1.1.lm1$fitted.values
#'     , GLM      =  Exam1.1.glm1$fitted.values
#'     , QB       =  Exam1.1.glm2$fitted.values
#'     , SM       =  Exam1.1.svyglm$fitted.values
#'   )
#' #-------------------------------------------------------------
#' ## One Method to plot  Figure1.1
#' #-------------------------------------------------------------
#' library(ggplot2)
#' 
#' Figure1.1   <-
#'   ggplot(
#'       data     = Newdata
#'     , mapping  = aes(x=x,y=y/Nx)
#'   )     +
#'   geom_point (
#'     mapping  = aes(colour="black")
#'   )  +
#'   geom_point (
#'     data     = Newdata
#'     , mapping  = aes(x=x,y=LM,colour="blue"),shape=2
#'   )  +
#'   geom_line(
#'     data     = Newdata
#'     , mapping  = aes(x=x,y=LM,colour="blue")
#'   )   +
#'   geom_point (
#'     data     = Newdata
#'     , mapping  = aes(x=x,y=GLM,colour="red"),shape=3
#'   ) +
#'   geom_smooth (
#'     data     = Newdata
#'     , mapping  = aes(x=x,y=GLM,colour="red")
#'     , stat     = "smooth"
#'   ) +
#'   theme_bw()    +
#'   scale_colour_manual (
#'     values=c("black","blue","red"),
#'     labels=c("observed","LM","GLM")
#'   )  +
#'   guides (
#'     colour   = guide_legend(title="Plot")
#'   ) +
#'   labs (
#'     title     = "Linear Model vs Logistic Model"
#'   ) +
#'   labs (
#'     y         = "p"
#'   )
#' print(Figure1.1)
#'
#'  #-------------------------------------------------------------
#' ## Another way to plot Figure 1.1
#' #-------------------------------------------------------------
#' newdata   <-
#'   data.frame(
#'     P     =  c(
#'                 Table1.1$y/Table1.1$Nx
#'               , Exam1.1.lm1$fitted.values
#'               , Exam1.1.glm1$fitted.values
#'                )
#'     , X     =  rep(Table1.1$x, 3)
#'     , group =  rep(c('Obs','LM','GLM'), each = length(Table1.1$x))
#'   )
#' 
#' Figure1.1      <-
#'   ggplot(
#'       data    = newdata
#'     , mapping = aes(x = X , y = P)
#'   )    +
#'   geom_point(
#'     mapping = aes(x = X , y = P, colour = group , shape=group)
#'   ) +
#'   geom_smooth(
#'     data    = subset(x = newdata, group == "LM")
#'     , mapping = aes(x=X,y=P)
#'     , col     = "green"
#'   ) +
#'   geom_smooth(
#'     data    = subset(x = newdata, group=="GLM")
#'     , mapping = aes(x = X , y = P)
#'     , col     = "red"
#'   ) +
#'   theme_bw() +
#'   labs(
#'     title   = "Linear Model vs Logistic Model"
#'   )
#' print(Figure1.1)
#'
#' #-------------------------------------------------------------
#' ## Correlation among p and fitted values using Gaussian link
#' #-------------------------------------------------------------
#' (lmCor <-
#'   cor(
#'     Table1.1$y/Table1.1$Nx,Exam1.1.lm1$fitted.values)
#' )
#' #-------------------------------------------------------------
#' ## Correlation among p and fitted values using quasi binomial link
#' #-------------------------------------------------------------
#' (glmCor  <-
#'   cor(
#'     Table1.1$y/Table1.1$Nx,Exam1.1.glm1$fitted.values)
#' )
NULL
