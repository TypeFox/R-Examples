predPow <-
function( data, obs = "observed", pred = "predicted", obs_mean = NULL, nu = 0, round = 4, extOut = FALSE, extOutFile = NULL ){
  mainfunc.stats( data, obs, pred, obs_mean, nu, round, "prediction", extOut, extOutFile, match.call() )
}

