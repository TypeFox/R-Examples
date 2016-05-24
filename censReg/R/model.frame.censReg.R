model.frame.censReg <- function( formula, ... ) {

   class( formula ) <- "lm"
   
   result <- model.frame( formula )

   return( result )
}

