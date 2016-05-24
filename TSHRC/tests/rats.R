
 library(TSHRC)
 data(rats)

 set.seed(42)
 twostage(rats$time, rats$delta, rats$group, nboot = 1e3)

