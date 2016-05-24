calibrationPower <-
function( data, obs = "observed", pred = "predicted", nu = 0, round = 4, extOut = FALSE, extOutFile = NULL ){
  mainfunc.stats( data, obs, pred, NULL, nu, round, "calibration", extOut, extOutFile, match.call() )
}

