## Wong (2010), Table 4.7 (p. 103), model 9

timings <- as.numeric(Sys.getenv("_R_CHECK_TIMINGS_"))
if(!is.na(timings) && timings > 60) {

library(logmult)
data(gss7590)

model <- rcL(gss7590, nd=2, weighting="none", se="jackknife", start=NA)

model
summary(model) # Jackknife standard errors are slightly different
               # from their asymptotic counterparts

stopifnot(all.equal(round(c(model$assoc$phi), d=3),
                    c(3.075, 3.474,  2.949, 2.460,
                      0.539, 1.686, -0.770, 1.891)))
stopifnot(all.equal(round(c(model$assoc$row[,,1]), d=3),
                    c(0.640,  0.239, -0.168, -0.711,
                      0.731, -0.217, -0.636,  0.121)))
stopifnot(all.equal(round(c(model$assoc$col[,,1]), d=3),
                    c(0.765,  0.250, -0.398, -0.273, -0.344,
                      0.071, -0.480, -0.198, -0.216,  0.824)))

# Check scores of heterogeneous model (scores not given by Wong)
model.heterog <- rcL(gss7590, nd=2, layer.effect="heterogeneous",
                     weighting="none", start=NA)
sep.models <- lapply(1:4, function(i) rc(gss7590[,,i], nd=2,
                                         weighting="none", start=NA))

# Level 2 does not converge properly because there are too few farmers:
# values are always slightly different and cannot be checked
stopifnot(isTRUE(all.equal(model.heterog$assoc$phi[-2,],
                           rbind(sep.models[[1]]$assoc$phi,
                                 sep.models[[3]]$assoc$phi,
                                 sep.models[[4]]$assoc$phi),
                           check.attributes=FALSE,
                           tolerance=1e-6)))
stopifnot(isTRUE(all.equal(model.heterog$assoc$row[,,-2],
                           array(c(sep.models[[1]]$assoc$row,
                                 sep.models[[3]]$assoc$row,
                                 sep.models[[4]]$assoc$row),
                                 dim=c(nrow(gss7590), 2, 3)),
                           check.attributes=FALSE,
                           tolerance=1e-6)))
stopifnot(isTRUE(all.equal(model.heterog$assoc$col[,,-2],
                           array(c(sep.models[[1]]$assoc$col,
                                 sep.models[[3]]$assoc$col,
                                 sep.models[[4]]$assoc$col),
                                 dim=c(ncol(gss7590), 2, 3)),
                           check.attributes=FALSE,
                           tolerance=1e-6)))

# Test anova
indep <- gnm(Freq ~ Education*Group + Occupation*Group, data=gss7590, family=poisson)
a <- anova(indep, model, model.heterog, test="LR")

} else {
cat("Skipped because _R_CHECK_TIMINGS_ is set and not empty.\n")
}
