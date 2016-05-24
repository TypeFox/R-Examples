fanSeries <- function(startValue, start, end, rates){
  ## Creates multivariate series that begin at start
  ## with the value startValue, and grow at the annual growth rates
  ## given by the vector rates until end. The resulting series will
  ## have length(rates) columns and the same frequency as start
  if(!inherits(start, "ti")) stop("start must be a ti")
  end <- ti(end, tif = tif(start))
  tvec <- seq(0, end - start - 1)/frequency(start)
  names(rates) <- paste(rates, "%", sep = "")
  tis(as.vector(startValue)*(1 + outer(tvec, rates/100)), 
     start = start)
}

tunnelSeries <- function(startValue,
                         start, 
                         end,
                         rate,
                         spreads){
  ## Creates a bivariate series with components that begin at
  ## start with the value startValue plus `spreads' percent, and
  ## grow at the annual growth rate given by rate until end. The
  ## resulting bivariate series will have the same frequency as start.
  ## spreads should be a vector of length 2.
  if(!inherits(start, "ti")) stop("start must be a ti")
  spreads <- rep(spreads, length.out = 2)
  offsets <- startValue*spreads/100
  midSeries <- fanSeries(startValue, start, end, rate)
  z <- cbind(midSeries + offsets[1], midSeries + offsets[2])
  dimnames(z) <- list(character(0),
					  paste(rate + spreads, "%", sep = ""))
  z
}
