library(bootES)
library(compiler)

## Constants
kFinal <- list(resultCols=c('stat', 'ci.low', 'ci.upper', 'bias', 'se'))

## Compiler Options
enableJIT(3)
. <- setCompilerOptions('optimize'=3)

## Functions to run the simulation
newResultMatrix <- function(reps) {
  cols <- c('stat', 'ci.low', 'ci.upper', 'bias', 'se')
  matrix(NA_real_, nrow=reps, ncol=length(cols), dimnames=list(NULL, cols))
}

doBootES <- function(R=1000, n, mean, sd) {
  dat <- data.frame(x=rnorm(n, mean, sd))
  res <- summary(bootES(dat, R, data.col="x"))
  as.vector(res)
}
doBootES <- cmpfun(doBootES)

## Parameters for the Monte Carlo Simulation
u1 <- 5
s1 <- 1
n <- 50

## Simulation 1
## Repetitions: 10,000
## Sample Size: 50
## bootES Samples: 1000
set.seed(1)
reps <- 1e4
dat <- data.frame(x=rnorm(n * reps, u1, s1))
mc.sim.1 <- sapply(seq_len(reps), function(x) {
  cat(x, "\n")
  idx <- (n * (x - 1) + 1):(n * (x - 1) + n)
  res <- summary(bootES(dat[idx, ,drop=FALSE], R=1000, data.col="x"))
  as.vector(res)
})
mc.sim.1 <- t(mc.sim.1)
colnames(mc.sim.1) <- kFinal$resultCols

cvg <- mc.sim.1[, 'ci.low'] < u1 & u1 < mc.sim.1[, 'ci.upper']
sum(cvg, na.rm=T)

## Simulation 2
## Repetitions: 10,000
## Sample Size: 50
## bootES Samples: 999
set.seed(2)
reps <- 1e4
mc.sim.2 <- newResultMatrix(reps)
for (i in seq_len(reps)) {
  cat(i, "\n")
  mc.sim.2[i, ] <- doBootES(R=999, n, u1, s1)
}
cvg <- mc.sim.2[, 'ci.low'] < u1 & u1 < mc.sim.2[, 'ci.upper']
sum(cvg, na.rm=T)

## Simulation 3
## Repetitions: 10,000
## Sample Size: 50
## bootES Samples: 1999
set.seed(2)
reps <- 1e4
mc.sim.3 <- newResultMatrix(reps)
for (i in seq_len(reps)) {
  cat(i, "\n")
  mc.sim.3[i, ] <- doBootES(R=1999, n, u1, s1)
}
cvg <- mc.sim.3[, 'ci.low'] < u1 & u1 < mc.sim.3[, 'ci.upper']
sum(cvg, na.rm=T)

## Simulation 3
## Repetitions: 10,000
## Sample Size: 50
## bootES Samples: 1999
set.seed(2)
u1 <- 5
s1 <- 1
reps <- 10000
mc.sim.4 <- newResultMatrix(reps)
for (i in seq_len(reps)) {
  cat(i, "\n")
  mc.sim.4[i, ] <- doBootES(R=1999, n=100, u1, s1)
}
cvg <- mc.sim.4[, 'ci.low'] < u1 & u1 < mc.sim.4[, 'ci.upper']
sum(cvg, na.rm=T)

mc.out <- "bootES.mc.RData"
if (!file.exists("bootES.mc.RData"))
  save(mc.sim.1, mc.sim.2, mc.sim.3, mc.sim.4, file="bootES.mc.RData",
       compress="xz")

## Formula for 'r' from 'd'
## r <- d / (sqrt(d^2 + 4))
## Simulation 4
## Repetitions: 10,000
## Sample Size: 50 x 50
## bootES Samples: 999

## Simulation Parameters
reps <- 1e4

## Population Parameters
# TODO
# u1 <- 5; s1 <- 1
# u2 <- 7; s2 <- 1
# cohens.d.sigma <- (u1 - u2) / s1
# r.pop <- cohens.d.sigma / sqrt(cohens.d.sigma^2 + 4)
# 
# set.seed(1)
# dat  <- data.frame(x=c(rnorm(50, u1, s1), rnorm(50, u2, s2)))
# dat$group <- rep(c('A', 'B'), each=50)
# mc.sim.4 <- newResultMatrix(reps)
# for (i in seq_len(reps)) {
#   cat(i, "\n")
#   mc.sim.4[i, ] <- doBootES(dat, R=999, 'x', 'group')
# }
# cvg <- mc.sim.4[, 'ci.low'] < u1 & u1 < mc.sim.4[, 'ci.upper']
# sum(cvg, na.rm=T)
# mc.sim.1 <- t(mc.sim.1)
# colnames(mc.sim.1) <- kFinal$resultCols
