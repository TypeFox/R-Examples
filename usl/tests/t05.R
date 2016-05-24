#
# Test function confint
#

library(usl)

data(raytracer)

options(digits=3, scipen=6)

set.seed(1103, kind = "default", normal.kind = "default")

u <- usl(throughput ~ processors, data = raytracer)

coef(u)

confint(u, parm=1)
confint(u, parm="sigma")

confint(u, parm=2)
confint(u, parm="kappa")

confint(u, parm=3)
confint(u, parm="none")
