.rebmix <- structure(list(
Preprocessing = c("histogram", "Parzen window", "k-nearest neighbour"),
Criterion = c("AIC", "AIC3", "AIC4", "AICc", "BIC", "CAIC", "HQC", "MDL2", "MDL5", "AWE", "CLC", "ICL", "ICL-BIC", "PC", "D", "SSE"),
Variables = c("continuous", "discrete"),
pdf = c("normal", "lognormal", "Weibull", "binomial", "Poisson", "Dirac", "gamma"),
pdf.nargs = c(2, 2, 2, 2, 1, 1, 2),
pdf.Variables = c("continuous", "continuous", "continuous", "discrete", "discrete", "discrete", "continuous"),
Restraints = c("rigid", "loose")),
.Names = c("Preprocessing", "Criterion", "Variables", "pdf", "pdf.nargs", "pdf.Variables", "Restraints"))

.rebmix.plot <- structure(list(
what = c("density", "marginal", "IC", "logL", "D", "distribution", "K")),
.Names = c("what"))

.rebmix.boot <- structure(list(
Bootstrap = c("parametric", "nonparametric")),
.Names = c("Bootstrap"))





