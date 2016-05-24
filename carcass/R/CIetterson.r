CIetterson <- function(s, s.lwr, s.upr, f, f.lwr, f.upr, J, s.time.variance="carcass age",
               f.time.variance="number of searches", nsim=1000, ci=0.95){
# either both s and f are dependent on calender date
# or both s is dependent on carcass age and f is dependent on the number of searches a carcass was exposed to

if(s.time.variance!="date") stopifnot(s.time.variance=="carcass age"&f.time.variance=="number of searches")
if(s.time.variance=="date") stopifnot(f.time.variance=="date")

s.a <- shapeparameter(s, s.lwr, s.upr)$a
s.b <- shapeparameter(s, s.lwr, s.upr)$b
f.a <- shapeparameter(f, f.lwr, f.upr)$a
f.b <- shapeparameter(f, f.lwr, f.upr)$b
n <- length(f)
N <- length(s)

p <- numeric(nsim)
for(r in 1:nsim) {
  sr <- rbeta(N, s.a, s.b)
  fr <- rbeta(n, f.a, f.b)
  if(N==1&n==1) p[r] <- ettersonEq14(s=sr, f=fr, J=J)
  if(N>1|n>1){
    if(s.time.variance!="date") p[r] <- ettersonEq14v2(s=sr, f=fr, J=J)
    if(s.time.variance=="date") p[r] <- ettersonEq14v1(s=sr, f=fr, J=J)
  }
  }
estp <- list(p.lower= quantile(p, prob=(1-ci)/2), p.upper=quantile(p, prob=1-(1-ci)/2))
return(estp)
}

