library(MonoPoly)

par <- 1
x <- seq(from=-1, to=1, length=21)
evalPol(x, par)

par <- rep(1, 2)
evalPol(x, par)

par <- 1:2
evalPol(x, par)

par <- 2:1
evalPol(x, par)

par <- rep(1, 3)
evalPol(x, par)

par <- 1:3
evalPol(x, par)

par <- 3:1
evalPol(x, par)

par <- rep(1, 4)
evalPol(x, par)

par <- 1:4
evalPol(x, par)

par <- 4:1
evalPol(x, par)
