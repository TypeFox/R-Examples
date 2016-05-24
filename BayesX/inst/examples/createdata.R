setwd("c:/arbeit/packages/BayesX")

# nonparametric function

x <- round(runif(300,-pi,pi),2)
time <- seq(1, 60)
time <- rep(time, 5)
y <- sin(x) + rnorm(300, 0, 0.3)
ytime <- sin(2*pi*time/60) + rnorm(300, 0, 0.3)

data <- data.frame(x,time,y,ytime)
write.table(data, "inst/examples/nonparametric.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)

# spatial effect
# surface
# sample path
# boundary file
# graph file
