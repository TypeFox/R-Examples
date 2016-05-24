
library("stsm")
library("KFKSDS")
library("numDeriv")

# derivative terms

xreg <- cbind(xreg = 1 * as.numeric(seq_len(length(JohnsonJohnson)) > 40))
m <- stsm::stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30, "xreg" = 3), 
  xreg = xreg)

#m@y[40] <- NA

# derivative terms

fnc.part <- function(x, model, type, i, j)
{
  m <- stsm::set.pars(model, x)
  m@y <- model@y - m@xreg %*% cbind(m@pars["xreg"])
  ss <- stsm::char2numeric(m, P0cov = FALSE)

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

g <- numDeriv::grad(func = fnc.part, x = m@pars, model = m, type = "v")
g
#[1]  0.2419351 -0.6267928  0.2972681 -2.7527843

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "apred", i = 1)
g
#[1]  -0.2571985   0.2439038  -0.1048057 -42.1365014

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "apred", i = 2)
g
#[1]  0.01526339  0.38288904 -0.19246241  0.88928568

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "aupd", i = 1)
g
#[1]  -0.2663266   0.2518944  -0.1081924 -43.1365014

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "aupd", i = 2)
g
#[1]  0.000680099 -0.227610738  0.113760420 -0.815062780

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "f", i = 2)
g
#[1] 162.7391 230.9178 190.3324   0.0000

g <- numDeriv::grad(func = fnc.part, x = m@pars, 
  model = m, type = "K", i = 1)
g
#[1] -0.2535469  0.6124994 -0.2893456  0.0000000

# KF.deriv

m2 <- m
m2@y <- m2@y - xreg * m2@pars["xreg"]
ss <- stsm::char2numeric(m2, P0cov = FALSE)
kf <- KF.deriv(m2@y, ss, xreg = m2@xreg)

colSums(kf$dv)
#       var1       var2       var3       xreg 
#  0.2419351 -0.6267928  0.2972681 -2.7527843 

g <- numDeriv::grad(func = fnc.part, x = m@pars, model = m, type = "v")
all.equal(g, as.vector(colSums(kf$dv)))
#[1] TRUE

colSums(kf$da.pred[,2,])
#        var1        var2        var3        xreg 
#  0.01526339  0.38288904 -0.19246241  0.88928568 

# KF.deriv.C

res <- KF.deriv.C(m2@y, char2numeric(m2, P0cov = FALSE), xreg = m@xreg)
colSums(res$dv)
#  0.2419351 -0.6267928  0.2972681 -2.7527843 
