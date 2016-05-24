#' Summary of DIF items
#' 
#' This function builds a table of DIF items specified in the \code{pwrrasch} object
#' 
#' @param object      \code{pwrrasch} object
#' @param all         If \code{TRUE}, all items are included in the table.
#' @param digits      Integer indicating the number of decimal places.
#' 
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # item parameters
#' ipar2 <- ipar1 <- seq(-3, 3, length.out = 20)
#' # model differential item function (DIF)
#' ipar2[10] <- ipar1[11]
#' ipar2[11] <- ipar1[10]
#' # simulation for b = 100 
#' simres <- pwr.rasch(100, ipar = list(ipar1, ipar2))
#' itemtable(simres)
#' }
itemtable <- function(object, all = FALSE, digits = 2) {
  
  if (length(object[[1]]) == 1) {
      
     ipar1 <- object$ipar[[1]]
     ipar2 <- object$ipar[[2]]
  
     min.items <- round(min(c(ipar1, ipar2)), digits = digits)
     max.items <- round(max(c(ipar1, ipar2)), digits = digits)
  
     n.DIF <- sum((ipar1 - ipar2) != 0)
  
  } else {
    
     object <- object[[1]]
    
     ipar1 <- object$ipar[[1]]
     ipar2 <- object$ipar[[2]]
    
     min.items <- round(min(c(ipar1, ipar2)), digits = digits)
     max.items <- round(max(c(ipar1, ipar2)), digits = digits)
    
     n.DIF <- sum((ipar1 - ipar2) != 0)
    
  }
  
  #------------------------------------------------------------------------------------------------------#
  
  if (all == FALSE) {
    
     pos <- which((ipar1 - ipar2) != 0)
    
     ipar1 <- ipar1[pos]
     ipar2 <- ipar2[pos]
    
  } else {
    
    pos <- 1:length(ipar1)
    
  }  
  
  if (length(pos) != 0) {
  
     itemtab <- data.frame(cbind(rbind(ipar1, ipar2, abs(ipar1 - ipar2)), 
                           c(sum(ipar1), sum(ipar2), sum(abs(ipar1 - ipar2)))))
           
     row.names(itemtab) <- c("ipar1", "ipar2", "DIF")
     names(itemtab) <- c(paste0("Pos", pos), "SUM")   
  
  } else {
    
    itemtab <- NULL
      
  }
  
  #------------------------------------------------------------------------------------------------------#
  # Output
  
  cat("\nSummary of DIF items specified in the pwrrasch object\n\n",
  
      "  item parameters: ", paste0("[", min.items, ", ", max.items, "]\n"),
      "  n items:         ", object$c, "\n",
      "  n DIF items:     ", n.DIF, "\n\n")
  
  if (length(pos) != 0) {
   
     printCoefmat(itemtab, digits = digits, has.Pvalue = FALSE)
     
     cat("\n")
    
  } else {
    
     cat("No DIF items specified.\n")
    
  }

  #------------------------------------------------------------------------------------------------------#
  
  return(invisible(itemtab))
  
}