#source("moninit")
library("RCALI")

# Cubature:
param <- list( output=1, input=2, dz=c(0,21), dp=c(100,0), tz=c(0,1))
califlopp("data2",c(fpollen,fseed),  param=param)
# Grid:
param <- list(method="grid",  grid=list(step=c(50,50)))
 califlopp("data2",dispf=c(3,1), param=param)
