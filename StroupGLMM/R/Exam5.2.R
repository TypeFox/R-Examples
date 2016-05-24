#' @title Example 5.2 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-164)
#' @name   Exam5.2
#' @docType data
#' @keywords datasets
#' @description Exam5.2 three factor main effects only design
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
#'    \code{\link{DataSet5.2}}
#'    
#' @importFrom lsmeans lsmeans contrast 
#' 
#' @examples
#' 
#'DataSet5.2$a <- factor( x = DataSet5.2$a)
#'DataSet5.2$b <- factor( x = DataSet5.2$b)
#'DataSet5.2$c <- factor(x  = DataSet5.2$c)
#' ##---first adding factor a in model
#' Exam5.2.lm1 <-
#'   lm(
#'       formula     = y~ a
#'     , data        = DataSet5.2
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
#' summary( Exam5.2.lm1 )
#'
#' library(lsmeans) 
#' ##---A first
#' ( Lsm5.2lm1    <-
#'   lsmeans::lsmeans(
#'      object  = Exam5.2.lm1
#'     , specs   = "a"
#'     # , ...
#'   )
#' )
#' ## lsmeans::contrast(object = Lsm5.2lm1 , method = "pairwise")
#' Anovalm1  <-   anova(object   = Exam5.2.lm1)
#' Anovalm1
#' 
#' ##---then adding factor b in model
#' Exam5.2.lm2 <-
#'   lm(
#'       formula     = y~ a + b
#'     , data        = DataSet5.2
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
#' summary( Exam5.2.lm1 )
#' (Lsm5.2lm2    <-
#'   lsmeans::lsmeans(
#'       object  = Exam5.2.lm2
#'     , specs   = "b"
#'     # , ...
#'   )
#' )
#' ## lsmeans::contrast(object = Lsm5.2lm2, method = "pairwise")
#' Anovalm2  <-   anova(object   = Exam5.2.lm2)
#' Anovalm2
#' 
#' ##---then adding factor c in model
#' Exam5.2.lm3 <-
#'   lm(
#'       formula     = y~ a + b + c
#'     , data        = DataSet5.2
#' #  , subset
#' #  , weights
#' #  , na.action
#'     , method      = "qr"
#'     , model       = TRUE
#' #  , x           = FALSE
#' #  , y           = FALSE
#'     , qr          = TRUE
#'     , singular.ok = TRUE
#'     , contrasts   = NULL
#' #  , offset
#' #  , ...
#'   )
#' summary( Exam5.2.lm3 )
#' (Lsm5.2lm3    <-
#'   lsmeans::lsmeans(
#'       object  = Exam5.2.lm3
#'     , specs   = "c"
#'     # , ...
#'   )
#' )
#' ## lsmeans::contrast(object = Lsm5.2lm3, method = "pairwise")
#' Anovalm3  <-  anova(object   = Exam5.2.lm3)
#' Anovalm3
#' 
#' ##---Now Change the order and add b first in model
#' Exam5.2.lm4 <-
#'   lm(
#'       formula     = y~  b
#'     , data        = DataSet5.2
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
#' summary( Exam5.2.lm4 )
#' (Lsm5.2lm4    <-
#'   lsmeans::lsmeans(
#'       object  = Exam5.2.lm4
#'     , specs   = "b"
#'     # , ...
#'   )
#' )
#' ## lsmeans::contrast(object = Lsm5.2lm4, method = "pairwise")
#' Anovalm4  <-  anova(object   = Exam5.2.lm4)
#' 
#' ##---then adding factor a in model
#' Exam5.2.lm5 <-
#'   lm(
#'       formula     = y~ b + a
#'     , data        = DataSet5.2
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
#' summary( Exam5.2.lm5 )
#' (Lsm5.2lm5    <-
#'   lsmeans::lsmeans(
#'       object  = Exam5.2.lm5
#'     , specs   = "a"
#'     # , ...
#'   )
#' )
#' ## lsmeans::contrast(object = Lsm5.2lm3, method = "pairwise")
#' Anovalm5  <-  anova(object   = Exam5.2.lm5)
#' Anovalm5
NULL
