ssim3p3 <-
function(truerate, n, r = 2, seed=NULL){
if (!is.null(seed)) {set.seed(seed)}
# n = number of simulations
# r = number of decimal places for the means
lp <- length(truerate)
# fmtd[i] <- freq dose i = mtd
fmtd <- rep(0, lp)
lastdose <- 0
# Matrix npt and ldt initialized to 0
mnpt <- mndlt <- matrix(rep(0, lp * n), lp, n)
# Progress bar creation
pb <- txtProgressBar(style=3)
setTxtProgressBar(pb,0)
for (i in 1:n){
sim <- sim3p3(truerate)
# Vector npt of the simulation i
mnpt[, i] <- sim$data$npt
# Vector dlt of the simulation i
mndlt[, i] <- sim$data$ndlt
if (sim$mtd %in% sim$data$dose) { 
lastdose <- lastdose + sim$lastdose
fmtd[sim$mtd] <- fmtd[sim$mtd] + 1
# Progress bar update
setTxtProgressBar(pb, i/n)
}
}
close(pb)
data <- CreData(lp)
data$npt <- round(apply(mnpt, 1, mean), r)
data$ndlt <- round(apply(mndlt, 1, mean), r)
pdlt <- round(apply(mndlt, 1, sum)/apply(mnpt, 1, sum), r)
exp = apply(mnpt, 1, sum) * 100 / sum(mnpt)
overshoot <- c(100, rep(0, lp))
for (i in 1:lp) {overshoot[i+1] <- overshoot[i] - exp [i]}
# pdlt[i] = probability to do a DLT at the dose i
list(data = cbind(data, pdlt, recommendation = fmtd * 100 / n, experimentation = round(exp, r), overshoot = round(overshoot[-1], r)), norecommendation = round(100 - sum(fmtd * 100 / n), r), mean.npt = sum(data$npt), mean.ndlt = sum(data$ndlt), mean.lastdose  = round(lastdose / n, r))
}
