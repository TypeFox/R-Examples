#' Birthweight data clustered on the mother. 
#'
#' This dataset is borrowed from "An introduction to Stata for health reserachers" (Juul and Frydenberg, 2010).
#' The dataset contains data on 189 mothers who have given birth to one or several children. In total, the dataset contains data on 487 births. 
#'
#' @docType data
#' @name clslowbwt
#' @usage data(clslowbwt)
#' @format The dataset is structured so that each row corresponds to one birth/child. It contains the following variables: 
#' \describe{
#'   \item{id}{the identification number of the mother.}
#'   \item{birth}{the number of the birth, i.e. "1" for the mother's first birth, "2" for the mother's second birth etc.}
#'   \item{smoke}{a categorical variable indicating if the mother is a smoker or not with levels "\code{0. No}" and "\code{1. Yes}".}
#'   \item{race}{the race of the mother with levels "\code{1. White}", "\code{2. Black}" or "\code{3. Other}".}
#'   \item{age}{the age of the mother at childbirth.}
#'   \item{lwt}{weight of the mother at last menstruational period (in pounds).}
#'   \item{bwt}{birthweight of the newborn.}
#'   \item{low}{a categorical variable indicating if the newborn is categorized as a low birthweight baby (<2500 grams) or not with levels "\code{0. No}" and "\code{1. Yes}".}
#'   \item{smoker}{a numeric indicator if the mother is a smoker or not. Recoded version of the variable "\code{smoke}" where "\code{0.No}" is recoded as "0" and "\code{1.Yes}" is recoded as "1".}
#'   \item{lbw}{a numeric indicator of whether the newborn is categorized as a low birthweight baby (<2500 grams) or not. Recoded version of the variable "\code{low}" where "\code{0.No}" is recoded as "0" and "\code{1.Yes}" is recoded as "1".}
#' }
#'
#'
#'\strong{The following changes have been made to the original data in Juul & Frydenberg (2010):}
#'
#' - The variable "\code{low}" is recoded into the numeric indicator variable "\code{lbw}": 
#' 
#' \code{clslowbwt$lbw <- as.numeric(clslowbwt$low == "1. Yes")}
#' 
#' - The variable "\code{smoke}" is recoded into the numeric indicator variable "\code{smoker}": 
#' 
#' \code{clslowbwt$smoker <- as.numeric(clslowbwt$smoke == "1. Yes")}
#' 
#' @references Juul, Svend & Frydenberg, Morten (2010). \emph{An introduction to Stata for health researchers}, Texas, Stata press, 2010 (Third edition).
#' @references \url{http://www.stata-press.com/data/ishr3.html}
#' 
             
NULL