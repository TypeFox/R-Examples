library(RTDE)


#####
# (1) small example

omega <- 1/2
m <- 10
n <- 100
obs <- cbind(rupareto(n), rupareto(n)) + rupareto(n)

#unit Pareto transform
z <- zvalueRTDE(obs, omega, nbpoint=m, output="relexcess")

MDPD(c(1/2, 1/4), dEPD, z$Z, alpha=0, rho=-1)
-sum(dEPD(z$Z, 1/2, 1/4, -1, log=TRUE))/m

MDPD(c(1/2, 1/4), dEPD, z$Z, alpha=0.05, rho=-1)
-21*sum(dEPD(z$Z, 1/2, 1/4, -1)^(0.05))/m + integrate(function(x) dEPD(x, 1/2, 1/4, -1)^(1.05), lower=1, upper=Inf)$value


#some check
dEPD(z$Z, 1/2, 1/4, -1)
do.call(dEPD, c(list(z$Z), as.list(c(1/2, 1/4)), rho=-1))
do.call(dEPD, c(list(z$Z), as.list(1/2), delta=1/4, rho=-1))
