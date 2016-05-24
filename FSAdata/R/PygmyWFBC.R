#' @title Biological data for Pygmy Whitefish from Dina Lake #1 (British Columbia), 2000 and 2001.
#' 
#' @description Biological information for Pygmy Whitefish (\emph{Prosopium coulterii}) from Dina Lake #1 (British Columbia), 2000 and 2001.
#' 
#' @name PygmyWFBC
#' 
#' @docType data
#' 
#' @format A data frame with 368 observations on the following 10 variables.
#'  \describe{
#'    \item{year}{Year of capture (2000, 2001)}
#'    \item{month}{Month of capture}
#'    \item{week}{Week within a month of capture}
#'    \item{net_no}{Unique net identification number}
#'    \item{fish_no}{Unique fish identification number}
#'    \item{fl}{Fork length (mm)}
#'    \item{wt}{Weight (g)}
#'    \item{sex}{Sex code (F=Female, M=Male, Imm=immature)}
#'    \item{mat}{Maturity code (Imm=immature, MG=maturing, MT=mature)}
#'    \item{scale_age}{Scale age (in years)}
#'    \item{oto_age}{Otolith age (in years)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Weight-Length
#'    \item Age Comparison
#'    \item Condition
#'    \item Length Frequency
#'  }
#'  
#' @concept 'Weight-Length' 'Age Comparison' Condition 'Length Frequency'
#' 
#' @source 2000 data from Appendix 9 of J. D. McPhail and R. J. Zemlak. October 2001. Pygmy whitefish studies on Dina Lake #1, 2000. Peace/Williston Fish and Wildlife Compensation Program, Report No. 245. 36pp plus appendices.  [Was (is?) from https://www.bchydro.com/pwcp/pdfs/reports/pwfwcp_report_no_245.pdf.]
#' 
#' 2001 data from Appendix 14 of Zemlak, R.J. and J.D. McPhail. 2003. Pygmy Whitefish studies on Dina Lake #1, 2001. Peace/Williston Fish and Wildlife Compensation Program Report No. 279. 35pp plus appendices.  [Was (is?) from http://www.bchydro.com/pwcp/pdfs/reports/pwfwcp_report_no_279.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(PygmyWFBC)
#' str(PygmyWFBC)
#' head(PygmyWFBC)
#' 
NULL