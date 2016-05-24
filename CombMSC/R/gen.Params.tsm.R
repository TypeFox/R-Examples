`gen.Params.tsm` <-
function(fmla, ...){
  make.Ind.Par <- function(size){ 2 * runif(size) - 1}
  list(AR.Pars = make.Ind.Par(as.numeric(fmla$model$AR)), MA.Pars = make.Ind.Par(as.numeric(fmla$model$MA)),
  SAR.Pars = make.Ind.Par(as.numeric(fmla$seasonal$model$AR)), SMA.Pars = make.Ind.Par(as.numeric(fmla$seasonal$model$MA)))
}

