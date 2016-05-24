#' @title Ages, lengths, and sexes of European Perch.
#' 
#' @description Assigned ages, measured fork lengths, and observed sexes for European Perch (\emph{Perca fluviatilis}) from Lake Tjuekemeer (The Netherlands).
#' 
#' @name EuroPerchTJ
#' 
#' @docType data
#' 
#' @format A data frame of 69 observations on the following 3 variables:
#'  \describe{
#'    \item{fl}{Fork lengths (cm).} 
#'    \item{age}{Assigned ages.}
#'    \item{sex}{Sex (female, male).} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Growth
#'    \item fon Bertalanffy
#'  }
#'  
#' @concept Growth 'von Bertalanffy'
#' 
#' @source From (approximately) Figure 2 in Mooij, W.M., J.M. Van Rooij, and S. Wijnhoven.  1999.  Analysis and comparison of fish growth from small samples of length-at-age data: Detection of sexual dimorphism in Eurasian perch as an example.  Transactions of the American Fisheries Society, 128:483-490.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(EuroPerchTJ)
#' str(EuroPerchTJ)
#' head(EuroPerchTJ)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(fl~age,data=EuroPerchTJ,subset=sex=="female",main="Female")
#' plot(fl~age,data=EuroPerchTJ,subset=sex=="male",main="Male")
#' par(op)
#' 
NULL