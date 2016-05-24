#
# test for the function 'mloglik.fd.deriv', it evaluates the analytical 
# expressions for the derivatives of the negative of the spectral 
# log-likelihood function 'mloglik.fd';
# numerical derivatives obtained with package 'numDeriv' are taken
# as benchmark
#
# NOTE: with 'cpar' and barrier term the analytical expression 
# and the numerical derivative do not match exactly
# use 'mcloglik.fd()' and 'mcloglik.fd.deriv()'

library("stsm")
library("numDeriv")

y <- log(AirPassengers)

pars <- c("var1" = 12, "var2" = 7, "var3" = 3.5, "var4" = 25)
nopars <- NULL

transP <- NULL #NULL #"StructTS" #"square"

bar <- list(type = "1", mu = 0)
#bar <- list(type = "1", mu = 0.01)
#bar <- list(type = "2", mu = 0.01)

m <- stsm.model(model = "BSM", y = y, pars = pars, nopars = nopars, 
  transPars = transP)

get.pars(m)
get.nopars(m)

mloglik.fd(model = m, barrier = bar)
mloglik.fd(x = m@pars, model = m, barrier = bar)

a1 <- mloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE, modcovgrad = TRUE,
  barrier = bar, version = "1")
a1$gradient
a1$hessian

a2 <- mloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE, modcovgrad = TRUE,
  barrier = bar, version = "2")
a2$gradient
a2$hessian

length(a1) == length(a2)
for (i in seq(along = a1))
  print(all.equal(a1[[i]], a2[[i]], check.attributes = FALSE))

g <- grad(func = mloglik.fd, x = m@pars, model = m, barrier = bar)
g
h <- hessian(func = mloglik.fd, x = m@pars, model = m, barrier = bar)
h

all.equal(a1$gradient, g, check.attributes = FALSE)
all.equal(a1$hessian, h, check.attributes = FALSE)

system.time(a1 <- mloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = FALSE, infomat = FALSE, modcovgrad = FALSE,
  barrier = bar, version = "1")$gradient)
system.time(a2 <- mloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = FALSE, infomat = FALSE, modcovgrad = FALSE,
  barrier = bar, version = "2")$gradient)
system.time(g <- grad(func = mloglik.fd, x = m@pars, model = m))
