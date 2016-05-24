
# tests for Kalman filter
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
convergence <- c(0.001, 10)

m@y[10] <- NA # NA before convergence
m@y[45] <- NA # NA after convergence

kf <- KF(m@y, ss)
kfc <- KF.C(m@y, ss, convergence = convergence)
kfd <- KF.deriv(m@y, ss, convergence = convergence)
kfdc <- KF.deriv.C(m@y, ss, convergence = convergence, return.all = TRUE)

# negative of the log-likelihood

kf$mll
kfc
kfd$mll
kfdc$mll

all.equal(kf$mll, kfc)
all.equal(kf$mll, kfd$mll)
all.equal(kf$mll, kfdc$mll)

#cbind(kfd$dv[,1], kfdc$dv[,1])

# derivative terms

fcn <- function(x, model, type, i, j)
{
  m <- stsm::set.pars(model, x)
  ss <- stsm::char2numeric(m)
  kf <- KF(m@y, ss)
  switch(type, 
    "v" = sum(kf$v, na.rm = TRUE), 
    "f" = sum(kf$f, na.rm = TRUE), 
    "K" = sum(kf$K[,i], na.rm = TRUE), 
    "apred" = sum(kf$a.pred[,i], na.rm = TRUE),
    "aupd" = sum(kf$a.upd[,i], na.rm = TRUE),
    "Ppred" = sum(kf$P.pred[i,j,], na.rm = TRUE),
    "Pupd" = sum(kf$P.upd[i,j,], na.rm = TRUE))
}

dv <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "v")
dv
colSums(kfd$dv, na.rm = TRUE)
colSums(kfdc$dv, na.rm = TRUE)
all.equal(dv, colSums(kfd$dv, na.rm = TRUE), check.attributes = FALSE)
all.equal(dv, colSums(kfdc$dv, na.rm = TRUE), check.attributes = FALSE)

df <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "f")
df
colSums(kfd$df, na.rm = TRUE)
colSums(kfdc$df, na.rm = TRUE)
all.equal(df, colSums(kfd$df, na.rm = TRUE), check.attributes = FALSE)
all.equal(df, colSums(kfdc$df, na.rm = TRUE), check.attributes = FALSE)

for (i in 1:4) {
dK <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "K", i = i)
print(all.equal(dK, colSums(kfd$dK[,i,]), check.attributes = FALSE))
print(all.equal(dK, colSums(kfdc$dK[,i,]), check.attributes = FALSE))
}

for (i in 1:4) {
daupd <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "aupd", i = i)
print(all.equal(daupd, colSums(kfd$da.upd[,i,]), check.attributes = FALSE))
#'da.upd' is not returned by KF.deriv.C
}

for (i in 1:4) {
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "apred", i = i)
print(all.equal(d, colSums(kfd$da.pred[,i,]), check.attributes = FALSE))
print(all.equal(d, colSums(kfdc$da.pred[,i,]), check.attributes = FALSE))
}

for (i in 1:4) for (j in 1:4) {
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "Ppred", i = i, j = j)
print(all.equal(d, colSums(kfd$dP.pred[i,j,,]), check.attributes = FALSE))
print(all.equal(d, colSums(kfdc$dP.pred[i,j,,]), check.attributes = FALSE))
}

for (i in 1:4) for (j in 1:4) {
d <- numDeriv::grad(func = fcn, x = m@pars, model = m, type = "Pupd", i = 1, j = 1)
print(all.equal(d, colSums(kfd$dP.upd[1,1,,]), check.attributes = FALSE))
#'dP.upd' is not returned by KF.deriv.C
}
