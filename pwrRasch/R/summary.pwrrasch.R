#' Object Summary
#' 
#' Generic \code{summary} function for the \code{pwrrasch} object
#' 
#' @param object      \code{pwrrasch} object
#' @param ...         Additional arguments affecting the summary produced.
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
#' ipar2[9] <- ipar1[12]
#' ipar2[12] <- ipar1[9]
#' # simulation for b = 100 
#' simres <- pwr.rasch(100, ipar = list(ipar1, ipar2))
#' summary(simres)
#' 
#' # item parameters
#' ipar2 <- ipar1 <- seq(-3, 3, length.out = 20)
#' # model differential item function (DIF)
#' ipar2[10] <- ipar1[11]
#' ipar2[11] <- ipar1[10]
#' # simulation for b = 100, 200, 300, 400, 500 
#' simres <- pwr.rasch(seq(100, 500, by = 100), ipar = list(ipar1, ipar2))
#' summary(simres)
#' }
summary.pwrrasch <- function(object, ...) {
  
  if (length(object[[1]]) == 1) {
     
     if ("type1" %in% names(object)) {
     
        with(object, {
    
        cat("\n Statistical Power Simulation for the Rasch model \n\n",
        
            "   b (number of persons in each group): ", b, "\n",
            "   c (number of items):                 ", c, "\n",            
            "   simulation runs:                     ", runs, "\n\n",
        
            "   Estimated statistical power: ", formatC(sum(H1.AC.p < sig.level) / runs, format = "f", digits = 3), "\n",
            "   Nominal significance level:  ", sig.level, "\n",
            "   Empirical significance level:", formatC(sum(H0.AC.p < sig.level) / runs, format = "f", digits = 3), "\n\n")
    
        })
            
     } else {
               
        with(object, {
                
          cat("\n Statistical Power Simulation for the Rasch model \n\n", 
        
              "   b (numer of persons in each group): ", b, "\n",
              "   c (numer of items):                 ", c, "\n",       
              "   simulation runs:                    ", runs, "\n\n",
        
              "   Estimated statistical power: ", formatC(sum(H1.AC.p < sig.level) / runs, format = "f", digits = 3), "\n",
              "   Nominal significance level:  ", sig.level, "\n\n")  
                  
        })
         
     }  

  #--------------------------------------------------------------------------------------------------------#   
  # length(object) > 1   
    
  } else {
    
     b <- unlist(lapply(object, function(x) x$b))
     c <- unique(unlist(lapply(object, function(x) x$c)))
     runs <- unique(unlist(lapply(object, function(x) x$runs)))
     pwr <- unlist(lapply(object, function(x) x$power))
     type1 <- unlist(lapply(object, function(x) x$type1))
    
     ###
    
     cat("\n Statistical Power Simulation for the Rasch model \n\n",
        
         "   b (numer of persons in each group): ", paste(b, collapse = ", "), "\n",
         "   c (numer of items):                 ", c, "\n",
         "   simulation runs:                    ", runs, "\n\n",
        
         "   Estimated statistical power: \n")
    
    for (i in 1:length(object)) {
            
      cat(paste0("      b = ", formatC(object[[i]]$b, digits = max(nchar(object[[i]]$b)) - 1, format = "d"), ":"),
          formatC(object[[i]]$power, format = "f", digits = 3) , "\n")
      
    }  
    
    ###
    
    if (!is.null(type1)) {
      
       cat("\n   Empirical significance level: \n")    
       for (i in 1:length(object)) {
        
           cat(paste0("      b = ", formatC(object[[i]]$b, digits = max(nchar(object[[i]]$b)) - 1, format = "d"), ":"),
               formatC(object[[i]]$type1, format = "f", digits = 3) , "\n")
        
       }  
      
    }
    
    cat("\n")      
    
  }      
     
}