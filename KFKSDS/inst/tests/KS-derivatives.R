
# tests for Kalman smoother
# compare numerical and analytical derivatives and 
# different funcions to compute the latter; 
# all of them should yield same results up to a tolerance

library("stsm")
library("KFKSDS")
library("numDeriv")

m <- stsm::stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
ss <- stsm::char2numeric(m)

#convergence <- c(0.001, length(m@y))
#convergence <- c(0.001, 10)

m@y[10] <- NA # NA before convergence
m@y[45] <- NA # NA after convergence

kf <- KF(m@y, ss, convergence = c(0.001, length(m@y)))
kf$convit
a <- KS(m@y, ss, kf)

kf <- KF(m@y, ss, convergence = c(0.001, 10))
kf$convit
b <- KS(m@y, ss, kf)

kfd <- KF.deriv(m@y, ss)
ksd <- KS.deriv(m@y, ss, kfd)

all.equal(a$r, b$r)
all.equal(a$N, b$N)
all.equal(a$ahat, b$ahat)
all.equal(a$varahat, b$varahat)
all.equal(a$r, ksd$r)
all.equal(a$N, ksd$N)
all.equal(a$ahat, ksd$ahat)
all.equal(a$varahat, ksd$varahat)

kf <- KF(m@y, ss, convergence = c(0.001, length(m@y)))
system.time(for(i in seq_len(10)) a <- KS(m@y, ss, kf))
kf <- KF(m@y, ss, convergence = c(0.001, 10))
system.time(for(i in seq_len(10)) b <- KS(m@y, ss, kf))

# derivative terms

fcn <- function(x, model, type, i, j)
{
  m <- stsm::set.pars(model, x)
  ss <- stsm::char2numeric(m)
  kf <- KF(m@y, ss)
  ks <- KS(m@y, ss, kf)
  switch(type, 
    "ahat" = sum(ks$ahat[,i], na.rm = TRUE), 
    "varahat" = sum(ks$varahat[i,j,], na.rm = TRUE),
    "r" = sum(ks$r[,i], na.rm = TRUE),
    "N" = sum(ks$N[i,j,], na.rm = TRUE))
}

for (i in 1:4)
{
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "ahat", i = i)
#colSums(ksd$dahat[,i,])
print(all.equal(d, colSums(ksd$dahat[,i,]), check.attributes = FALSE))
}

for (i in 1:4)
{
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "r", i = i)
print(all.equal(d, colSums(ksd$dr[,i,]), check.attributes = FALSE))
}

for (i in 1:4) for (j in 1:4)
{
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "varahat", i = i, j = j)
print(all.equal(d, colSums(ksd$dvarahat[,i,j,]), check.attributes = FALSE))
}

for (i in 1:4) for (j in 1:4)
{
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "N", i = i, j = j)
print(all.equal(d, colSums(ksd$dN[,i,j,]), check.attributes = FALSE))
}
