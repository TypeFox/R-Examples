cat("\ntest HLfit:")
# HLfit

data(wafers)
## Gamma GLMM with log link
hl <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
            resid.model = ~ X3+I(X3^2) ,data=wafers)

expect_equal(hl$APHLs$p_v,-1157.609,tolerance=1e-03)
