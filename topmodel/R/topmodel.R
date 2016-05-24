topmodel <- function(parameters, topidx, delay, rain, ET0, verbose = F, Qobs = NA) {

  ## deal with verbosity:

  v = 1
  if(verbose && is.na(Qobs)) v <- 6

  ## number of iterations

  if (is.vector(parameters)) iterations <- 1
  else {
    if (is.matrix(parameters)) iterations <- dim(parameters)[1]
    else return(NA)
  }

  ## check data inconsistencies

  if(length(parameters)/iterations < 11) stop("Incorrect number of parameters")
  if(length(rain[is.na(rain)]) != 0 ) stop("Rain contains incorrect values such as NA's")

  ## check dangerous parameter values (not implemented)

  ## check whether the function should return E or for Qsim
  ## Adjust lengthResult accordingly

  if(length(Qobs) == 1 && is.na(Qobs)) {
    Qobs <- -9999    # go for Qsim
    lengthResult <- length(rain)*iterations
  }
  else {
    if( length(Qobs[Qobs>=0])!=length(Qobs) || length(Qobs)!= length(rain)) {
      print("Check Qobs for negative values or wrong length")
      return(NA)
    }
    else {
      Qobs[is.na(Qobs)] <- -1
      lengthResult <- iterations
    }
  }

  ## running the model...

  result <- .C("topmodel",
               PACKAGE = "topmodel",
               as.double(t(parameters)),
               as.double(as.matrix(topidx)),
               as.double(as.matrix(delay)),
               as.double(rain),
               as.double(ET0),
               as.double(Qobs),
               as.integer(length(as.double(as.matrix(topidx)))/2),
               as.integer(length(rain)),
               as.integer(iterations),
               as.integer(length(delay[,1])),
               as.integer(v),
               result = double(v * lengthResult))$result

  ## formatting of the results

  if(v == 6) {
    result <- matrix(result,ncol=6)
    result <- list(
                   Q  = matrix(result[,1], ncol=iterations),
                   qo = matrix(result[,2], ncol=iterations),
                   qs = matrix(result[,3], ncol=iterations),
                   S  = matrix(result[,4], ncol=iterations),
                   fex= matrix(result[,5], ncol=iterations),
                   Ea = matrix(result[,6], ncol=iterations)
                   )
  }

  if((Qobs == -9999) && (iterations > 1) && (v == 1)) result <- matrix(result, ncol= iterations)

  return(result)

}
