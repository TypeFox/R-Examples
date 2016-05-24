# a logistic transition function for simulationg data from a STCC-GARCH
   tr.func   <- function(tr.par, tr.var){
      g <- tr.par[1]
      cons <- tr.par[2]
      (1.0 + exp(-g*(tr.var-cons)))^(-1)
   }
