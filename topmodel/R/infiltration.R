infiltration <- function(rain, parameters) {

  result <- .C("infiltration",
               PACKAGE = "topmodel",
               as.double(rain),
               as.double(parameters),
               as.integer(length(rain)),
               result = double(length(rain)))$result

  result[result < -9000] <- NA
  
  return(result)

}
