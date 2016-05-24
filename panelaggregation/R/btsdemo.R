#' Randomly Generated Panel Dataset
#'
#' This data was created by simulation to mimmick a firm level dataset stemming from business tendency surveys.
#' The data was simulated because of privacy concerns with micro level firm data. For convenience the dataset 
#' contains two different date notations. Also 5 qualitative 3-item questions are included. Business tendency survey data
#' is often weighted with company size represented by the number of employees. Thus the weight column is quantitative and
#' its distribution is somewhat (!) reasonable with respect to the distribution of employees in a typical firm sample. 
#'
#'
#' \itemize{
#'   \item uid unique company identifier
#'   \item year numeric year column
#'   \item weight quantitative weight
#'   \item question\_1
#'   \item question\_2
#'   \item question\_3
#'   \item question\_4
#'   \item question\_5
#'   \item group group to mimmick different sectors / branches of trade
#'   \item altGroup another alternative grouping columns
#'   \item sClass a column denoting discrete size classes small (S), medium (M) and large (L)
#'   \item date\_qtrly quarterly dates stored in a single column. 
#' }
#'
#' @author Matthias Bannert
#' @format A data frame with 27000 rows and 13 variables
#' @source Randomly generated in R using the sample generator from https://github.com/mbannert/gateveys/blob/master/R/gateveys.R
#' @name btsdemo
NULL