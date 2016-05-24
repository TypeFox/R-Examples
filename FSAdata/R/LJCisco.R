#' @title Ages and lengths of Longjaw Cisco from two locations in Lake Michigan.
#' 
#' @description Assigned age (by scales) and total length of Longjaw Cisco (\emph{Leucichthys alpenae}) captured at two locations in Lake Michigan.
#' 
#' @name LJCisco
#' 
#' @docType data
#' 
#' @format A data frame with 378 observations on the following 3 variables.
#'  \describe{
#'    \item{age}{Assigned age (by scales).}
#'    \item{tl}{Measured total length (mm).} 
#'    \item{loc}{Capture location (\code{NE}=northeast and \code{S}=south).} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item von Bertalanffy
#'  }
#'  
#' @concept Growth 'von Bertalanffy'
#' 
#' @source Simulated from age-length data provided in tables 2 and 3 of Jobes, F.W.  1946.  The age, growth, and distribution of the longjaw cisco, \emph{Leucichthys alpenae} Koelz, in Lake Michigan.  Transactions of the American Fisheries Society.  76:215-247.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LJCisco)
#' str(LJCisco)
#' head(LJCisco)
#' op <- par(mfrow=c(1,2))
#' plot(tl~age,data=LJCisco,subset=loc=="NE",main="northeast")
#' plot(tl~age,data=LJCisco,subset=loc=="S",main="south")
#' par(op)
#' 
NULL