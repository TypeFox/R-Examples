#' @title Example1.2 from Generalized Linear Mixed Models: Modern Concepts, Methods and Applications by Walter W. Stroup(p-9)
#' @name   Exam1.2
#' @docType data
#' @keywords datasets
#' @description Exam1.2 is used to see types of model effects by plotting regression data
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
#'    \code{\link{Table1.2}}

#' @importFrom ggplot2 ggplot
#' 
#' @examples
#' 
#' #-------------------------------------------------------------
#' ## Plot of multi-batch regression data discussed in Article 1.3
#' #-------------------------------------------------------------
#' data(Table1.1)
#' Table1.2$Batch <- factor(x  = Table1.2$Batch)
#' library(ggplot2)
#' Plot  <-
#'  ggplot(
#'    data    = Table1.2
#'    , mapping = aes(y = Y, x =X,colour=Batch,shape=Batch)
#'  )      +
#'  geom_point() +
#'  geom_smooth(
#'    method  = "lm"
#'    , fill    =  NA
#'  ) +
#'  labs(
#'    title   = "Plot of Multi Batch Regression data"
#'  )        +
#'  theme_bw()
NULL
