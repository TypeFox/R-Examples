simCrm <-
function(prior, firstdose = NA, truerate = prior, cohortsize = 3, target = 1/3, nptmax = 24, nmaxmtd = nptmax, nmaxdose = nptmax, sd = 1.34, approach = "bayes", model = "power", method = "fpost", nextlevel = "ntarget", upskipping = F, downskipping = F, lastdose = NA, graphic = F, seed=NULL){#browser()
if (!is.null(seed)) {set.seed(seed)}
# prob = probability vector
data <- CreData(length(prior))
if (graphic == T)
mprob <- matrix(prior, length(prior))
dlt <- dose <- vector()
if (firstdose %in% data$dose)
nextdose <- firstdose
else if (firstdose %in% NA)
nextdose <- Crm(data, prior, target, nptmax, nmaxmtd, nmaxdose, sd, approach, model, method, nextlevel)$nextdose
else
stop("first specified not available.")
while (nextdose %in% data$dose) {
lastdose <- nextdose
dose <- c(dose, lastdose)
dlt <- c(dlt, rbinom(1, cohortsize, truerate[lastdose]))
data <- updata(data, lastdose, cohortsize, dlt[length(dlt)])
sim <- Crm(data, prior, target, nptmax, nmaxmtd, nmaxdose, sd, approach, model, method, nextlevel, upskipping, downskipping, lastdose)
nextdose <- sim$nextdose
if (graphic == T)
mprob <- cbind(mprob,sim$prob)
}
prob <- sim$prob
if (graphic == T) {
palette(gray(seq(0,.9,len=ncol(mprob))))
colors <- palette(gray(seq(0,.9,len=ncol(mprob))))
nc <- length(colors)
plot(data$dose,prior,type="b",pch=15,xlab="dose level",ylab="DLT probability",col=colors[nc],lwd=2,ylim=c(0,1))
abline(h=target,col="gray",lty=3)
for (i in 1:nc - 1) {
lines(data$dose,mprob[,i + 1],type="b",col=colors[nc - i],pch=15,lwd=2)
}
}
mtd <- sim$mtd
list(data = data, dose = dose, ndlt = dlt, mtd = mtd, lastdose = lastdose, prob = prob)
}
