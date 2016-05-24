`fit.Models.tsm` <-
function(fmla, data.Vector, ...){
  arima(data.Vector, order = c(fmla$model$AR, fmla$model$I, fmla$model$MA),
  seasonal = list(order=c(fmla$seasonal$model$AR,  fmla$seasonal$model$I,
  fmla$seasonal$model$MA), period = fmla$seasonal$period))
}

