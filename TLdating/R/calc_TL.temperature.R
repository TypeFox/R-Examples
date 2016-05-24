#' calculate temperature vector
#'
#' This function calculates the temperature vector.
#'
#' @param nPoints
#'  \link{numeric} (\bold{required}): number of points.
#' @param Tmax
#'  \link{numeric} (\bold{required}): Maximum temperature.
#' @param Hrate
#'  \link{numeric} (\bold{required}): Heating rate.
#' @param an_temp
#'  \link{numeric} (with default): Annealing temperature.
#' @param an_time
#'  \link{numeric} (with default): Annealing time.
#' @param rec_ramp2PH
#'  \link{logical} (with default): Indicate if the signal was record during the ramp up to the preheat temperature.
#' @param rec_duringPH
#'  \link{logical} (with default): Indicate if the signal was record during the preheat plateau.
#'
#'
#' @return
#'  This function provides a new \code{\linkS4class{TLum.Results}} object containing the times and temperature vectors.
#'
#'
#' @author David Strebler, University of Cologne (Germany).
#'
## @export calc_TL.temperature

calc_TL.temperature <- function(

  nPoints,
  Tmax,
  Hrate,
  an_temp=0,
  an_time=0,
  rec_ramp2PH = FALSE,
  rec_duringPH= FALSE

  ){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(nPoints)){
    stop("[calc_TL.temperature] Error: Input 'nPoints' is missing.")
  }else if (!is.numeric(nPoints)){
    stop("[calc_TL.temperature] Error: Input 'nPoints' is not of type 'numeric'.")
  }
  if (missing(Tmax)){
    stop("[calc_TL.temperature] Error: Input 'Tmax' is missing.")
  }else if (!is.numeric(Tmax)){
    stop("[calc_TL.temperature] Error: Input 'Tmax' is not of type 'numeric'.")
  }
  if (missing(Hrate)){
    stop("[calc_TL.temperature] Error: Input 'Hrate' is missing.")
  }else if (!is.numeric(Hrate)){
    stop("[calc_TL.temperature] Error: Input 'Hrate' is not of type 'numeric'.")
  }

  if (!is.numeric(an_temp)){
    stop("[calc_TL.temperature] Error: Input 'an_temp' is not of type 'numeric'.")
  }
  if (!is.numeric(an_time)){
    stop("[calc_TL.temperature] Error: Input 'an_time' is not of type 'numeric'.")
  }

  if(!is.logical(rec_ramp2PH) || is.na(rec_ramp2PH)){
    stop("[calc_TL.temperature] Error: Input 'rec_ramp2PH' is not of type 'logical'.")
  }

  if(!is.logical(rec_duringPH) || is.na(rec_duringPH)){
    stop("[calc_TL.temperature] Error: Input 'rec_duringPH' is not of type 'logical'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check


  if(an_temp<0){
    an_temp <- 0
  }else if(an_temp>Tmax){
    an_temp <- Tmax
  }

  # ------------------------------------------------------------------------------

  Dtot <- an_time + Tmax/Hrate
  Dramp1 <- an_temp/Hrate
  Dplateau <- an_time
  Dramp2 <- (Tmax-an_temp)/Hrate

  if(an_temp==0){
    D1 <- 0
    D2 <- 0
    D3 <- Tmax/Hrate

  } else {

    if(rec_ramp2PH){
      D1 <- Dramp1
    }else{
      D1 <- 0
    }

    if(rec_duringPH){
      D2 <- an_time
    }else{
      D2 <- 0
    }

    D3 <- Dramp2
  }

  Drec <- D1+D2+D3
  Dstep <- Drec/nPoints
  Tstep <- Hrate*Dstep

  full.times <- seq(from=Dstep, to=Dtot, by=Dstep)
  rec.times <- seq(from=Dstep, to=Drec, by=Dstep)

  full.temperatures <- vector()
  rec.temperatures <- vector()

  full.temperatures[1] <- Tstep

  if(rec_ramp2PH){
    rec.temperatures[1] <- Tstep
  }else{
    rec.temperatures[1] <- an_temp+Tstep
  }

  j <- 1

  for(i in 2: length(full.times)){
    if(full.times[i]< Dramp1){
      full.temperatures[i] <- full.temperatures[i-1]+Tstep

      if(rec_ramp2PH){
        j <- j+1
        rec.temperatures[j] <- full.temperatures[i]
      }

    }else if(full.times[i] < Dramp1 + Dplateau){
      full.temperatures[i] <- an_temp

      if(rec_ramp2PH){
        j <- j+1
        rec.temperatures[j] <- full.temperatures[i]
      }

    }else if(full.times[i] < Dtot){
      full.temperatures[i] <- full.temperatures[i-1]+Tstep

      j <- j+1
      rec.temperatures[j] <- full.temperatures[i]
    }else{
      full.temperatures[i] <- Tmax

      j <- j+1
      rec.temperatures[j] <- full.temperatures[i]
    }
  }

  result <- list(temperatures=rec.temperatures,
                 times = rec.times,
                 full.temperatures=full.temperatures,
                 full.times =full.times
                 )

  new.TLum.Results.calc_TL.temperature <- set_TLum.Results(data = result)

  return(new.TLum.Results.calc_TL.temperature)
}
