## Wong (2010), Table 2.7 (p. 48-49)
# See also ?rc

timings <- as.numeric(Sys.getenv("_R_CHECK_TIMINGS_"))
if(!is.na(timings) && timings > 60) {

library(logmult)
data(gss8590)

# The table used in Wong (2010) is not perfectly consistent
# with that of Wong (2001)
tab <- margin.table(gss8590[,,c(2,4)], 1:2)
tab[2,4] <- 49

model <- rc(tab, nd=2, weighting="none", se="jackknife", start=NA)

model
summary(model) # Jackknife standard errors are slightly different
               # from their asymptotic counterparts

# Compare with bootstrap standard errors
model2 <- rc(tab, nd=2, weighting="none", se="bootstrap", start=NA)
plot(model, conf.ellipses=0.95)
summary(model2)

# A few scores differ from reported results by .001
stopifnot(isTRUE(all.equal(round(c(model$assoc$phi), d=3),
                           c(2.601, 1.522))))
stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                           c(0.743,  0.088, -0.200, -0.632,
                            0.276, -0.061, -0.777,  0.562))))
stopifnot(isTRUE(all.equal(round(c(model$assoc$col), d=3),
                           c(0.765,  0.020, -0.322, -0.55, 0.088,
                            -0.137, -0.549, -0.120, -0.01, 0.816))))

} else {
cat("Skipped because _R_CHECK_TIMINGS_ is set and not empty.\n")
}
