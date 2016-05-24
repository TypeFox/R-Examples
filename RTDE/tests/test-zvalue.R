library(RTDE)


#####
# (1) small example

omega <- 1/2
m <- 10
n <- 100
obs <- cbind(rupareto(n), rupareto(n)) + rupareto(n)

#unit Pareto transform
zvalueRTDE(obs, omega, output="orig")

relexcess(zvalueRTDE(obs, omega, output="orig"), m)
zvalueRTDE(obs, omega, nbpoint=m, output="relexcess")

#unit Frechet transform
zvalueRTDE(obs, omega, output="orig", marg="ufrechet")

relexcess(zvalueRTDE(obs, omega, output="orig", marg="ufrechet"), m)
zvalueRTDE(obs, omega, nbpoint=m, output="relexcess", marg="ufrechet")

