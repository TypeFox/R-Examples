`gen.Data.tsm` <-
function(fmla, pars, data.Size, rand.gen, test.Size, ...){
  dat <- sarima.Sim(round(data.Size/fmla$seasonal$period), period = fmla$seasonal$period, model = list(ar = pars$AR.Pars,
  ma = pars$MA.Pars, order = c(fmla$model$AR, fmla$model$I, fmla$model$MA)),
  seasonal = list(ar = pars$SAR.Pars, ma = pars$SMA.Pars, order = c(fmla$seasonal$model$AR,
  fmla$seasonal$model$I, fmla$seasonal$model$MA), period=fmla$seasonal$period))
  temp <- splitTrainTest(dat, numTrain = length(dat) - test.Size)
  list(data.Vector = temp$train, test.Vector = temp$test)
}

