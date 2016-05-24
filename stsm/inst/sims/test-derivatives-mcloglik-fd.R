#
# test for the function 'mloglik.fd.deriv', it evaluates the analytical 
# expressions for the derivatives of the negative of the spectral 
# log-likelihood function 'mloglik.fd';
# numerical derivatives obtained with package 'numDeriv' are taken
# as benchmark
#
# NOTE: the analytical derivatives for 'mcloglik.fd()' 
# are not currently implemented for a model parameterized 
# according to a non-null 'transPars'
# a barrier term is not considered either

library("stsm")
library("numDeriv")

y <- log(AirPassengers)

pars <- c("var2" = 7, "var3" = 3.5, "var4" = 25)
cpar <- c("var1" = 12)
nopars <- NULL

transP <- NULL

bar <- list(type = "1", mu = 0)

m <- stsm.model(model = "BSM", y = y, 
  pars = pars, cpar = cpar, nopars = nopars, 
  transPars = transP)

get.pars(m)
get.cpar(m)
get.nopars(m)

mcloglik.fd(model = m, barrier = bar)
mcloglik.fd(x = m@pars, model = m, barrier = bar)

a1 <- mcloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE)
a1$gradient
a1$hessian

a2 <- mcloglik.fd.deriv(model = m, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE)
a2$gradient
a2$hessian

length(a1) == length(a2)
for (i in seq(along = a1))
  print(all.equal(a1[[i]], a2[[i]], check.attributes = FALSE))

g <- grad(func = mcloglik.fd, x = m@pars, model = m, barrier = bar)
g
h <- hessian(func = mcloglik.fd, x = m@pars, model = m, barrier = bar)
h

all.equal(a1$gradient, g, check.attributes = FALSE)
all.equal(a1$hessian, h, check.attributes = FALSE)
