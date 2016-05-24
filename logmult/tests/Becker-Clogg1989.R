## Becker & Clogg (1989), Table 5 (p. 145)
# See also ?rc

timings <- as.numeric(Sys.getenv("_R_CHECK_TIMINGS_"))
if(!is.na(timings) && timings > 60) {

library(logmult)
data(color)

# "Uniform weights" in the authors' terms mean "no weighting" for us,
# and "average marginals" means "marginal" with rcL
# See ?rc for "marginals"
caithness.unweighted <- rc(color[,,1], nd=2, weighting="none",
                           se="jackknife", start=NA)
caithness.marginal <- rc(color[,,1], nd=2, weighting="marginal",
                         se="jackknife", start=NA)
aberdeen.unweighted <- rc(color[,,2], nd=2, weighting="none",
                          se="jackknife", start=NA)
aberdeen.marginal <- rc(color[,,2], nd=2, weighting="marginal",
                        se="jackknife", start=NA)

caithness.unweighted
caithness.marginal
aberdeen.unweighted
aberdeen.marginal

se(caithness.unweighted)
se(caithness.marginal)
se(aberdeen.unweighted)
se(aberdeen.marginal)

summary(caithness.unweighted)
summary(caithness.marginal)
summary(aberdeen.unweighted)
summary(aberdeen.marginal)

stopifnot(isTRUE(all.equal(c(round(caithness.unweighted$assoc$phi, d=3)),
                           c(3.067, 0.375))))
stopifnot(isTRUE(all.equal(c(round(caithness.unweighted$assoc$row, d=3)),
                           c(0.411,  0.478, -0.122, -0.767, 0.472,
                            -0.030, -0.804,  0.362))))
stopifnot(isTRUE(all.equal(c(round(caithness.unweighted$assoc$col, d=3)),
                           c(0.574, 0.279,  0.121, -0.261, -0.714,
                             0.271, 0.148, -0.838,  0.449, -0.030))))

stopifnot(isTRUE(all.equal(c(round(caithness.marginal$assoc$phi, d=3)),
                           c(0.494, 0.110))))
stopifnot(isTRUE(all.equal(c(round(caithness.marginal$assoc$row, d=3)),
                           c(0.894,  1.045, -0.166, -1.52, 1.246,
                             0.271, -1.356,  0.824))))
stopifnot(isTRUE(all.equal(c(round(caithness.marginal$assoc$col, d=3)),
                           c(c(1.313, 0.425, -0.019, -1.213, -2.567,
                               0.782, 0.518, -1.219,  0.945,  0.038)))))

# This is one of the very few articles reporting jackknife standard errors
# Our results are very close for the unweighted solution;
# but for marginal weighting their errors are probably wrong,
# as they are larger than the unweighted ones, while phi are much smaller.
stopifnot(isTRUE(all.equal(c(round(se(caithness.unweighted)$phi, d=3)),
                           c(0.282, 0.053))))
stopifnot(isTRUE(all.equal(c(round(se(aberdeen.unweighted)$phi, d=3)),
                           c(0.221, 0.040))))



## Same with rc(M)-L model
# See also ?rcL
data(color)

# "Uniform weights" in the authors' terms mean "no weighting" for us,
# and "average marginals" means "marginal" with rcL
# See ?rc for "marginals"
unweighted <- rcL(color, nd=2, weighting="none",
                  layer.effect="heterogeneous", se="jackknife", start=NA)
marginal <- rcL(color, nd=2, weighting="marginal",
                layer.effect="heterogeneous", se="jackknife", start=NA)
unweighted
marginal

# (our standard errors are much smaller for the marginal-weighted case)
summary(unweighted)
summary(marginal)

opar <- par(mfrow=c(1, 2))
plot(marginal, layer="Caithness", conf.ellipses=0.95)
plot(marginal, layer="Aberdeen", conf.ellipses=0.95)
par(opar)

stopifnot(isTRUE(all.equal(c(round(unweighted$assoc$phi, d=3)),
                           c(3.067,  2.854, 0.375, 0.294))))
stopifnot(isTRUE(all.equal(c(round(unweighted$assoc$row, d=3)),
                           c(0.411, 0.478, -0.122, -0.767, 0.472,
                            -0.030, -0.804, 0.362, 0.492, 0.399, -0.128,
                            -0.763, 0.452, -0.053, -0.797, 0.397))))
stopifnot(isTRUE(all.equal(c(round(unweighted$assoc$col, d=3)),
                           c(0.574, 0.279,  0.121, -0.261, -0.714, 0.271, 0.148, 
                            -0.838, 0.449, -0.030,  0.595,  0.225, 0.135, -0.231,
                            -0.724, 0.420,  0.094, -0.867,  0.206, 0.147))))

stopifnot(isTRUE(all.equal(c(round(marginal$assoc$phi, d=3)),
                           c(0.465, 0.415, 0.111, 0.085))))
stopifnot(isTRUE(all.equal(c(round(marginal$assoc$row, d=3)),
                           c(0.883, 1.035, -0.188, -1.554,  1.244,  0.268,
                            -1.337, 0.863,  1.175,  0.936, -0.246, -1.512,
                             1.184, 0.184, -1.317, 0.976))))
stopifnot(isTRUE(all.equal(c(round(marginal$assoc$col, d=3)),
                           c(1.361, 0.424, -0.051, -1.301, -2.731,  0.803,  0.544,
                            -1.175, 0.970,  0.075,  1.398,  0.212, -0.067, -1.260,
                            -2.848, 0.811,  0.463, -1.177,  0.942, 1.144))))

# This is one of the very few articles reporting jackknife standard errors
# Our results are very close for the unweighted solution;
# but for marginal weighting their errors are probably wrong,
# as they are larger than the unweighted ones, while phi are much smaller.
stopifnot(isTRUE(all.equal(c(round(se(unweighted)$phi, d=3)),
                           c(c(0.282, 0.221, 0.053, 0.04)))))
} else {
cat("Skipped because _R_CHECK_TIMINGS_ is set and not empty.\n")
}
