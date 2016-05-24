cat("\ntest of inverse Gamma:")

data(wafers)
HLig <- HLfit( y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),
       family=Gamma(log),rand.family=inverse.Gamma(log),
       resid.model= ~ X3+I(X3^2) ,data=wafers)
expect_equal(HLig$APHLs$p_v,-1157.523,tolerance=1e-3)