#' @title Catch and effort for three Snapper species in a depletion experiment.
#' 
#' @description Catch and effort for three Snapper species (\emph{Pristipomoides zonatus}, \emph{Pristipomoides auricilla}, and \emph{Etelis carbunculUs}) in a depletion experiment around Pathfinder Reef in the Mariana Archipelago.
#' 
#' @name Pathfinder
#' 
#' @docType data
#' 
#' @format A data frame with 13 observations on the following 5 variables.
#'  \describe{
#'    \item{date}{Date (1984)}
#'    \item{effort}{Fishing effort (line-hours of a bottom hand-line)}
#'    \item{Pzonatus}{Catch of \emph{Pristipomoides zonatus}} 
#'    \item{Pauricilla}{Catch of \emph{Pristipomoides auricilla}}
#'    \item{Ecarbunculus}{Catch of \emph{Etelis carbunculUs}}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{ 
#'    \item Depletion methods 
#'    \item Leslie method 
#'    \item DeLury method 
#'    \item Population size 
#'    \item Abundance 
#'    \item Catchability 
#'  }
#'  
#' @concept Abundance 'Population Size' Leslie DeLury Depletion Catchability
#' 
#' @source From Table 3 of Polovina, J.J.  1985.  A variable catchability version of the Leslie model with application to an intensive fishing experiment on a multispecies stock.  Fishery Bulletin 84:423-428.  [Was (is?) from https://swfsc.noaa.gov/publications/CR/1986/8679.PDF.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Pathfinder)
#' str(Pathfinder)
#' head(Pathfinder)
NULL