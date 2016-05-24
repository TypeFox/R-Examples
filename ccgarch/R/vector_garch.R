#computing vector GARCH volatilities: valid for CCC, ECCC, DCC, EDCC models

   vector.garch <- function(dvar, a, A, B){
      dvar <- dvar^2           # dvar = eps
     .Call("vector_garch", dvar, a, A, B)
   }
