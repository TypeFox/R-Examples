"Add2" <-
function(nx, ny, X, Y, quantile, alternative) 
{

pxI <- (X+0.5)/(nx+1)
pyI <- (Y+0.5)/(ny+1)

nxI <- nx+1
nyI <- ny+1

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

