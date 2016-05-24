#' @title Evaluation of the Orientation Services 2002
#' @description A sample survey on students evaluation of the Orientation services was conducted across the
#' 13 Faculties of University of Naples - Federico II in five waves: participants were asked to express their ratings
#'  on a 7 point scale (1 = "very unsatisfied", 7 = "extremely satisfied").
#' Here only the dataset collected during 2002 is loaded.
#' @aliases univer
#' @usage data(univer)
#' @format The description of subjects' covariates is:
#' \describe{
#' \item{\code{Faculty}}{A factor variables, with levels ranging from 1 to 13 indicating the coding 
#' for the different university faculties}
#' \item{\code{Freqserv}}{A factor with levels: 0 = for not regular users, 1 = for regular users}
#' \item{\code{Age}}{Variable indicating the age of the respondent in years}
#' \item{\code{Gender}}{A factor with levels: 0 = man, 1 = woman}
#' \item{\code{Diploma}}{A factor with levels:  1 = classic studies, 2 = scientific studies, 3 = linguistic,
#'  4 = Professional, 	5 = Technical/Accountancy, 6 = others}
#' \item{\code{Residence}}{A factor with levels: 1 = city NA, 2 = district NA, 3 = others}
#'  \item{\code{ChangeFa}}{A factor with levels: 1 = changed faculty, 2 = not changed faculty}
#'	}
#' Analyzed ordinal variables (Likert ordinal scale):
#'  1 = "extremely unsatisfied", 2 = "very unsatisfied", 3="unsatisfied", 4="indifferent", 5="satisfied", 6 = "very satisfied",
#' 7="extremely satisfied"
#' \describe{
#' \item{\code{Informat}}{Level of satisfaction about the collected information}
#' \item{\code{Willingn}}{Level of satisfaction about the willingness of the staff}
#' \item{\code{Officeho}}{Judgment about the Office hours}
#' \item{\code{Competen}}{Judgement about the competence of the staff}  
#' \item{\code{Global}}{Global satisfaction}
#' }
#' @keywords datasets
#' @details 
#' \describe{
#' Period of data collection: 2002 \cr
#' Mode of collection: questionnaire \cr
#' Number of observations: 2179 \cr
#' Number of subjects' covariates: 7 \cr
#' Number of analyzed items: 5
#' }
#' @source \url{http://www.labstat.it/home/wp-content/uploads/2015/09/univer02.txt}


"univer"

