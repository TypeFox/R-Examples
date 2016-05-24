#' @title Glycogen content in rat livers
#'
#' @description An experiment measuring the glycogen content in rat
#' livers that uses subsampling.
#'
#' @details Six rats were randomised to three treatment conditions
#' (two per condition). Their livers were divided into three pieces
#' and two measurements were taken on each piece for 6 x 3 x 2 = 36
#' observations. The rats are the experimental units and there are two
#' levels of subsampling.
#'
#' @format A data frame with 36 rows and 4 variables:
#' \describe{
#' 
#'   \item{Glycogen:}{Glycogen content.}
#' 
#'   \item{Treatment:}{A numeric treatment group indicator (1, 2, 3).}
#' 
#'   \item{Rat:}{Rat identification number. The numbers are not
#' unique; each treatment has a rat numbered 1 and 2, but these not
#' the same rats (rats are nested under treatment).}
#' 
#'   \item{Liver:}{A numeric variable indicating the liver piece. Each
#' liver was divided into three pieces.}  }
#'
#' @references Sokal RR and Rohlf FJ (1995). \emph{Biometry}. WH
#' Freeman and Co., New York, NY.
#' 
#' @docType data
#' @name glycogen
NULL
