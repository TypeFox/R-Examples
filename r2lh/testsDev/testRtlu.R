cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                      Begin testRtlu                     +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/rtlu.R")

r2lUnivFactor(f2,"toto")
r2lUnivOrdered(o2,"toto")
r2lUnivDiscrete(d2,"toto")
r2lUnivContinuous(c2,"toto")


rtlu(f1)
rtlu(f2)

rtlu(o1)
rtlu(o2)

rtlu(i1)
rtlu(i2)

rtlu(n1)
rtlu(n2)


cat("\n---------------------------------------------------------------
---                        End testRtlu                     ---
---------------------------------------------------------------\n")

