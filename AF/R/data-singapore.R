#' Case-control study on oesophageal cancer in Chinese Singapore men. 
#'
#' This dataset is borrowed from "Aetiological factors in oesophageal cancer in Singapore Chinese" by De Jong UW, Breslow N, Hong JG, Sridharan M, Shanmugaratnam K (1974).
#'
#' @docType data
#' @name singapore
#' @usage data(singapore)
#' @format The dataset contains the following variables:
#' \describe{
#' \item{Age}{age of the patient.} 
#' \item{Dial}{dialect group where 1 represent "\code{Hokhien/Teochew}" and 0 represent "\code{Cantonese/Other}".}
#' \item{Samsu}{a numeric indicator of whether the patient consumes Samsu wine or not.}
#' \item{Cigs}{number of cigarettes smoked per day.}
#' \item{Bev}{number of beverage at "burning hot" temperatures ranging between 0 to 3 different drinks per day.}
#' \item{Everhotbev}{a numeric indicator of whether the patients ever drinks "burning hot beverage" or not. Recoded from the variable "\code{Bev}".}
#' \item{Set}{matched set identification number.}
#' \item{CC}{a numeric variable where 1 represent if the patient is a case, 2 represent if the patient is a control from the same ward as the case and 3 represent if the patient is control from orthopedic hospital.}
#' \item{Oesophagealcancer}{a numeric indicator variable of whether the patient is a case of oesophageal cancer or not.} 
#' }
#' 
#' \strong{The following changes have been made to the data from the original data in De Jong UW (1974):}
#' 
#' - The variable "\code{Bev}" is recoded into the numeric indicator variable "\code{Everhotbev}":
#' 
#' \code{singapore$Everhotbev <- ifelse(singapore$Bev >= 1, 1, 0)}
#' 
#' @references 	De Jong UW, Breslow N, Hong JG, Sridharan M, Shanmugaratnam K. (1974). Aetiological factors in oesophageal cancer in Singapore Chinese. \emph{Int J Cancer} Mar 15;13(3), 291-303.
#' @references \url{http://faculty.washington.edu/heagerty/Courses/b513/WEB2002/datasets.html}

NULL
