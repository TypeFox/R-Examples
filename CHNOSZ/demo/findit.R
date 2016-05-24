## findit() calculations for sulfur species
basis("CHNOS+")
basis("pH", 5)
species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", "S3O6-2",
  "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
# to minimize the standard deviations of the 
# logarithms of activity the species
objective <- "SD"
# the variables we are interested in
vars <- list(O2=c(-50, -15), pH=c(0, 14), T=c(275, 375))
# optimize logfO2 at constant T and pH
f1 <- findit(vars[1], objective, T=325, P=350, niter=3)
title("S.D. of equilibrium log activities of sulfur species")
# optimize logfO2 and pH at constant T
f2 <- findit(vars[1:2], objective, T=325, P=350, res=16, niter=5)
title("S.D. of equilibrium log activities of sulfur species")
# optimize logfO2, pH and T (at constant P ...)
f3 <- findit(vars, objective, P=350, res=10, niter=10)
title("S.D. of equilibrium log activities of sulfur species")
# the results
print(f1.out <- sapply(f1$value, tail, 1))
print(f2.out <- sapply(f2$value, tail, 1))
print(f3.out <- sapply(f3$value, tail, 1))
# with more variables, we should find a greater degree of optimization
stopifnot(f2.out["SD"] < f1.out["SD"])
stopifnot(f3.out["SD"] < f2.out["SD"])
