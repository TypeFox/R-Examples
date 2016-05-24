"Add4" <-
function(nx, ny, X, Y, quantile, alternative) 
{

pxI <- (X+1)/(nx+2)
pyI <- (Y+1)/(ny+2)

nxI <- nx+2
nyI <- ny+2
estI <- pxI - pyI
stderr <- sqrt( pxI*(1-pxI)/nxI + pyI*(1-pyI)/nyI )

if(alternative=="two.sided")
 {
  lower <- estI - quantile*stderr
  upper <- estI + quantile*stderr
 }

if(alternative=="less")
 {
  lower <- (-1)
  upper <- estI + quantile*stderr
 }

if(alternative=="greater")
 {
  lower <- estI + quantile*stderr
  upper <- 1
 }

px <- X/nx
py <- Y/ny
estimate <- px-py

list(conf.int=c(lower,upper),
estimate=estimate)
   
}

