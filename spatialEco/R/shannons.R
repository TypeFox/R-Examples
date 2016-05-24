#' @title Shannon's Diversity (Entropy) Index
#' @description Calculates Shannon's Diversity Index and Shannon's Evenness Index
#'
#' @param x        data.frame object containing species counts
#' @param counts   Are data counts (TRUE) or relative proportions (FALSE) 
#' @param margin   Calculate diversity for rows or columns. c("row", "col")    
#'
#' @return data.frame with columns "H" (Shannon's diversity) and "evenness" (Shannon's evenness where H / max( sum(x) ) ) 
#'
#' @note The sdist argument will produce an evenly spaced sample, whereas n produces a fixed sized sample. The p (proportional) argument calculates the percent of the line-length. 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references
#' Shannon, C. E. and W. Weaver (1948) A mathematical theory of communication. The Bell System Technical Journal, 27:379-423.
#' Simpson, E. H. (1949) Measurement of diversity. Nature 163:688
#' Roth, D. S., I. Perfecto, and B. Rathcke (1994) The effects of management systems on ground-foraging ant diversity in Costa Rica. Ecological Applications 4(3):423-436.
#'
#' @examples
#' # Using Costa Rican ant diversity data from Roth et al. (1994)
#' data(ants)
#'   
#' # Calculate diversity for each covertype ("col") 
#' shannons(ants[,2:ncol(ants)], counts = FALSE, margin = "col")
#'
#' # Calculate diversity for each species ("row") 
#' ant.div <- shannons(ants[,2:ncol(ants)], counts = FALSE, margin = "row")
#'   names(ant.div) <- ants[,1]
#'   ant.div
#'	
#' @export
shannons <- function(x, counts = TRUE, margin = "row") {
  if( margin == "row") { m = 1 } else { m = 2 }
  s <- function(x) { -1 * sum( x[!is.na(x) > 0] * log(x[!is.na(x) > 0]) ) } 
  x.sum <- apply(x, MARGIN = m, sum, na.rm=TRUE)
  if( counts ) { x <- sweep(x, MARGIN = m, x.sum, "/") }
  H <- apply(x, MARGIN = m, FUN = s)
    if(counts) {
      EH <- H / log(max(x.sum))
      return( data.frame(H=H, evenness=EH) )
	  } else {
	  return( H )  
    }	
}
