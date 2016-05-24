#' @title Catch and effort data for Deepwater Caridean Shrimp.
#' 
#' @description Catch and effort data for Deepwater Caridean Shrimp (\emph{Heterocarpus laevigatus}) from 15 days in 1984 from near Alamagan Islan in the Marian Archipelago (near Guam).
#' 
#' @details Catch (kg) and effort (trap-nights) of Deepwater Caridean Shrimp (\emph{Heterocarpus laevigatus}) from 15 days in 1984 from near Alamagan Islan in the Marian Archipelago (near Guam).  The data start on 9-Jan-1984.  Catches were recorded separately for standard traps and in pyramid traps.
#' 
#' The original authors estiamted populations size using the Leslie method with the cumulative catch from the combined catch in the standard and pyramid traps, but with a CPE computed from just the catch in standard traps.
#' 
#' @name ShrimpGuam
#' 
#' @docType data
#' 
#' @format A data frame with 15 observations on the following 4 variables.
#'  \describe{
#'    \item{day}{Day of the catch.  Day 9 is 9-Jan-1984.} 
#'    \item{standard}{Catch (kg) of of Shrimp in the standard traps.} 
#'    \item{pyramid}{Catch (kg) of of Shrimp in the pyramid traps.} 
#'    \item{effort}{Total effort (trap-nights) for the standard traps.} 
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
#' @source From Table 1 of Ralson, S. 1986.  An intensive fishing experiment for the Caridean Shrimp, \emph{Heterocarpus laevigatus}, at Alamagn Island in the Mariana Archipelago.  Fishery Bulletin 84:927-934.  [Was (is?) from http://fishbull.noaa.gov/844/ralston.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(ShrimpGuam)
#' str(ShrimpGuam)
#' head(ShrimpGuam)
#' 
#' ## Computations by the original authors
#' # CPE for just the standard traps
#' ShrimpGuam$CPE <- ShrimpGuam$standard/ShrimpGuam$effort
#' # Total catch in both traps
#' ShrimpGuam$total <- ShrimpGuam$standard+ShrimpGuam$pyramid
#' # Cumulative catch in both traps (with the Ricker modification)
#' ShrimpGuam$cumCatch <- cumsum(ShrimpGuam$total)-ShrimpGuam$total/2
#' 
NULL