Xsc.onesample <-
function(data, pi0){
if(missing(data) || missing(pi0))
stop("data and/or pi0 missing.")

nreads <- rowSums(data)
fit.MoM <- DM.MoM(data)
MoM.theta <- fit.MoM$theta
MoM.pi <- fit.MoM$pi

rank.Bj <- sum(pi0>0) - 1
Xsc <- Xsc.statistics(MoM.pi, MoM.theta, nreads, pi0)[1]
p.value <- 1-pchisq(q=Xsc, df=rank.Bj, ncp=0, lower.tail=TRUE)

RAD.mean.test <- list(Xsc, p.value)
names(RAD.mean.test) <- c("Xsc statistics", "p value")

return(RAD.mean.test)
}
