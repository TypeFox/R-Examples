heckit <- function( selection, outcome,
                    data = sys.frame(sys.parent()), method="2step", ... ) {

   thisCall <- match.call()
   result <- selection( selection = selection, outcome = outcome,
      data = data, method = method, ... )

   result$call <- thisCall

   return( result )
}
