"Wald" <-
function(nx, ny, X, Y, quantile, alternative )
{

px <- X/nx
py <- Y/ny
estimate <- px-py

stderr<- sqrt( px*(1-px)/nx  +  py*(1-py)/ny )

if(alternative=="two.sided")
 {
  lower<-estimate - quantile*stderr
  upper<-estimate + quantile*stderr
 }

if(alternative=="less")
 {
  lower<- (-1)
  upper<-estimate + quantile*stderr
 }

if(alternative=="greater")
 {
  lower<-estimate + quantile*stderr
  upper<- 1
 }

list(conf.int=c(lower,upper),
estimate=estimate)
   
}

