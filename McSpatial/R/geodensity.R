geodensity <- function(longvar,latvar,window=.25,kern="tcub",alldata=FALSE) {

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,.5,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { .75*(1-psi^2) } }
  if (kern=="bisq")  { wgt <- function(psi) { (15/16)*((1-psi^2)^2) } }
  if (kern=="tcub")  { wgt <- function(psi) { (70/81)*((1 - abs(psi)^3)^3) } }
  if (kern=="trwt")  { wgt <- function(psi) { (35/32)*((1 - psi^2)^3) } }

  longvar <- 2*pi*longvar/360
  latvar  <- 2*pi*latvar/360

  if (alldata==FALSE) {
    fit <- locfit(~lp(longvar,latvar,nn=window,deg=1)) 
    xev <- lfeval(fit)$xev
    nt = length(xev)/2
    target <- t(array(xev,dim=c(2,nt)))
  }
  if (alldata==TRUE) {
    target <- cbind(longvar,latvar)
    nt = nrow(target)
  }

# pmin function assures that dist is either <1 or truly equal to 1
  dens.target <- array(0,dim=nt)
  for (i in seq(1:nt)) {
    dist <- pmin(sin(latvar)*sin(target[i,2]) + cos(latvar)*cos(target[i,2])*cos(target[i,1]-longvar),  1)
    dist <- acos(dist)*3958
    h = quantile(dist,window)
    k <- ifelse(dist<h, wgt(dist/h), 0) 
    dens.target[i] = mean(k)/h
  }

  if (alldata==FALSE) {
    denshat <- smooth12(target,dens.target,cbind(longvar,latvar))
  }
  if (alldata==TRUE) {denshat <- dens.target}

  target <- 360*target/(2*pi)
  out <- list(target,dens.target,denshat) 
  names(out) <- c("target","dens.target","denshat")
  return(out)    
}

