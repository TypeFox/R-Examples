"cf.forecasts" <-
function(m1,m2)
  { rmse <- rmse(m1,m2)
    mae <- mae(m1,m2)
    return(c(rmse,mae)) }

