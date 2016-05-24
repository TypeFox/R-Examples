### Detect the tobit model type based on how selection and
### outcome formulas are set up.
### returns the model type

detectModelType <- function(selection, outcome) {
   type <- 0
   if(!inherits( selection, "formula" )) {
      stop( "argument 'selection' must be a formula in function 'selection()'" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula in function 'selection()'" )
   }
                                        #
   if(inherits(outcome, "formula")) {
      if( length( outcome ) != 3 ) {
         stop( "argument 'outcome' must be a 2-sided formula in function 'selection()'" )
      }
      type <- 2
   }
   else if(inherits(outcome, "list")) {
      if(length(outcome) == 1) {
         outcome <- outcome[[1]]
         type <- 2
      }
      else if(length(outcome) == 2) {
         if(inherits(outcome[[1]], "formula")) {
            if( length( outcome[[1]] ) != 3 ) {
               stop( "argument 'outcome[[1]]' must be a 2-sided formula in function 'selection()'" )
            }
         }
         else
             stop( "argument 'outcome[[1]]' must be a formula in function 'selection()'" )
         if(inherits(outcome[[2]], "formula")) {
            if( length( outcome[[2]] ) != 3 ) {
               stop("argument 'outcome[[2]]' must be a 2-sided formula in function 'selection()'" )
            }
         }
         else
             stop( "argument 'outcome[[2]]' must be a formula in function 'selection()'" )
         type <- 5
      }
      else
          stop("argument 'outcome' must contain 1 or 2 components in function 'selection()'")
   }
   else
       stop("argument 'outcome' must be either a formula or a list of two formulas in function 'selection()'" )
   return(type)
}
