`sarima.Sim` <-
function (n = 20, period = 12, model, seasonal, rand.Gen.Fun = rnorm, rand.Gen.Seas = rnorm)
{
  tmp <- arima.sim( model, n*period , rand.gen = rand.Gen.Fun)
  tmp <- tmp[!is.na(tmp) & tmp != 0]
  for (i in 1:period) {
    tmp2 <- arima.sim( seasonal, n , rand.gen = rand.Gen.Seas)
    tmp2 <- tmp2[!is.na(tmp2) & tmp2 != 0]
    tmp[i + period * 0:(n-1)] <-
      tmp[i + period * 0:(n-1)] + tmp2
  }
  ts(tmp, frequency=period)
}

