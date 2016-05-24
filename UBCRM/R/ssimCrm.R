ssimCrm <-
function(prior, n, firstdose = NA, truerate = prior, cohortsize = 3, target = 1/3, nptmax = 24, nmaxmtd = nptmax, nmaxdose = nptmax, sd = 1.34, approach = "bayes", method = "fpost", model = "power", nextlevel = "ntarget", upskipping = F, downskipping = F, r = 2, seed=NULL){
# N = number of simulations
if (!is.null(seed)) {set.seed(seed)}
lp <- length(prior)
fmtd <- rep(0, lp)
# fmtd[i] <- freq dose k = mtd
lastdose <- vector()
# Matrix npt and ldt initialized to 0
mnpt <- mndlt <- mprob <- matrix(rep(0, lp * n), lp, n)
# Progress bar creation
pb <- txtProgressBar(style=3)
setTxtProgressBar(pb,0)
for (i in 1:n){
sim <- simCrm(prior, firstdose, truerate, cohortsize, target, nptmax, nmaxmtd, nmaxdose, sd, approach, model, method, nextlevel, upskipping, downskipping, lastdose)
# Vector npt of the simulation i
mnpt[, i] <- sim$data$npt
# Vector dlt of the simulation i
mndlt[, i] <- sim$data$ndlt
# Vector probability of the simulation i
mprob[, i] <- sim$prob
if (sim$mtd %in% sim$data$dose) { 
lastdose <- c(lastdose, sim$lastdose)
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
list(data = cbind(data, pdlt, recommendation = fmtd * 100 / n, experimentation = round(exp, r), overshoot = round(overshoot[-1], r)), norecommendation = 100 - sum(fmtd * 100 / n), mean.npt = sum(data$npt), mean.ndlt = sum(data$ndlt), mean.lastdose  = round(mean(lastdose), r), mean.prob = apply(mprob, 1, mean))
}
