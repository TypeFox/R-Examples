working.smoothvariogram <-
function(coords,gamma.long,dd,nmin=3)
{
  ## coords    : coordinates missing bins
  ## gamma : variogram in long form
  ## dd    : distance and time lag difference
  ## nmin  : number of minimal neighbors used for interpolation
  ## value : the weighted averages of the variogram 
  coord.s <- (gamma.long$s)
  coord.t <- (gamma.long$t)
  out <- numeric()
  for(i in 1:nrow(coords)){
	   tmps <- abs(coord.s -coords[i,1])/dd[1]
	   tmpt <- abs(coord.t - coords[i,2])/dd[2]
	   tmp <- pmax(tmps,tmpt)
       ii <- order(tmp)[1:nmin]
	   np <- gamma.long[ii,]$n
	   g <- gamma.long[ii,]$gamma
       tmpdd <- sqrt( tmps[ii]^2+tmpt[ii]^2) 
       ww <- np/tmpdd
       ww <- ww / sum(ww)
       out <- rbind(out,sum(ww*g))
  }
  out
}
