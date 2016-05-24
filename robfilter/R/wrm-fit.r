WRMfit <- function(xdat, ydat, x0=xdat[length(xdat)], weight.vec=rep(1,length(xdat))){
  weight.vec     <- weight.vec/sum(weight.vec)
  xdat        <- xdat[weight.vec!=0]
  ydat        <- ydat[weight.vec!=0]
  weight.vec     <- weight.vec[weight.vec!=0]
  sm          <- rep(0,length(xdat))
  for (i in 1:length(xdat)){
     sm[i] <- (weighted.quantile.top((ydat[i]-ydat[-i])/(xdat[i]-xdat[-i]), weight.vec[-i], 0.5) + weighted.quantile((ydat[i]-ydat[-i])/(xdat[i]-xdat[-i]),weight.vec[-i],0.5))/2
  }
  bm <- (weighted.quantile.top(sm, weight.vec, 0.5) + weighted.quantile(sm, weight.vec, 0.5))/2
  am <- (weighted.quantile.top(ydat-bm*(xdat-x0), weight.vec, 0.5) + weighted.quantile(ydat-bm*(xdat-x0), weight.vec, 0.5))/2
  c(am,bm)  # am is fitted value
}
