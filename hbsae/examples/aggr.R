d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop)

# by default aggregate over all areas
global <- aggr(sae)
EST(global); SE(global)

# aggregation to broad area
# first build aggregation matrix
M <- d$Xpop[, c("area22", "area23", "area24")] / d$Xpop[, "(Intercept)"]
M <- cbind(1 - rowSums(M), M); colnames(M)[1] <- "area21"
est.area2 <- aggr(sae, M)
EST(est.area2); SE(est.area2)
COV(est.area2)  # covariance matrix
