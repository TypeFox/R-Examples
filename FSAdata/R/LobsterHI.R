#' @title Catch and effort data for Hawaiian Islands Slipper Lobster.
#' 
#' @description Catches of Slipper Lobster (\emph{Scyllarides squammosus}) in three categories from the vicinity of Laysan Bank, Hawaiian Islands on 34 consecutive days in 1986.
#' 
#' @details Catch (numbers) of lobster in three categories - legal (tail weight greater than 85g), sublegal (tail weight less than 85g), and berried (egg-bearing females).  Sublegal and berried lobsters were returned to the water.
#' 
#' The vessel fished between 11 June and 14 July 1986 in the vicinity of Laysan Island and its associated bank.  The daily operations of the vessel involved deploying and hauling 1,125 Fathom Plus lobster traps set in strings spaced at 30 m intervals.  They were fished in 7 strings of about 160 traps each and baited with Pacific Mackerel, \emph{Scomber japonicus}. Strings were soaked overnight and retrieved the following day; therefore, the standard unit of effort is the trap-haul.
#' 
#' @name LobsterHI
#' 
#' @docType data
#' 
#' @format A data frame with 34 observations on the following 6 variables.
#'  \describe{
#'    \item{day}{Day of the catch} 
#'    \item{legal}{Number of legal lobsters caught.} 
#'    \item{sublegal}{Number of sub-legal lobsters caught.}
#'    \item{berried}{Number of egg-bearing lobsters caught.} 
#'    \item{total}{Total number of lobsters caught.} 
#'    \item{effort}{Total daily effort expended.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population size 
#'    \item Abundance 
#'    \item Depletion methods
#'    \item Leslie method
#'    \item DeLury method 
#'    \item Catchability
#'  }
#'  
#' @concept Abundance 'Population Size' Leslie DeLury Depletion Catchability
#' 
#' @source From Table 1 of Clarke, R.P., and S.S. Yoshimoto.  1990.  Application of the Leslie model to commercial catch and effort of the slipper lobster, \emph{Scyllarides squammosus}, fishery in the northwestern hawaiian islands.  Marine Fisheries Review, 52(2):1-7.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LobsterHI)
#' str(LobsterHI)
#' head(LobsterHI)
#' 
NULL