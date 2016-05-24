cleanrf <- function(mpcross)
{
  if (is.null(mpcross$rf) | (sum(is.na(mpcross$rf$theta))==0)) 
	return(colnames(mpcross$finals))

  rf <- mpcross$rf$theta
  nm <- apply(rf, 2, function(x) sum(is.na(x)))

  while (sum(is.na(rf))>0)
  {
    m <- which.max(nm)
    rf <- rf[-m, -m]
    nm <- apply(rf, 2, function(x) sum(is.na(x)))
  }

  return(colnames(rf))
}
