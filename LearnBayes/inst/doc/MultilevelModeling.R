### R code from vignette source 'MultilevelModeling.Rnw'

###################################################
### code chunk number 1: MultilevelModeling.Rnw:17-24
###################################################
d <- data.frame(Name=c("Clemente", "Robinson", "Howard", "Johnstone",
      "Berry",  "Spencer",  "Kessinger", "Alvarado", "Santo", 
      "Swaboda", "Petrocelli",  "Rodriguez", "Scott",  "Unser", 
      "Williams",  "Campaneris",  "Munson",   "Alvis"),
      Hits=c(18, 17, 16, 15, 14, 14, 13, 12, 11, 
           11, 10, 10, 10, 10, 10,  9,  8,  7),
      At.Bats=45)


###################################################
### code chunk number 2: MultilevelModeling.Rnw:59-64
###################################################
library(LearnBayes)
laplace.fit <- laplace(betabinexch, 
                       c(0, 0),
                       d[, c("Hits", "At.Bats")])
laplace.fit


###################################################
### code chunk number 3: MultilevelModeling.Rnw:68-73
###################################################
mcmc.fit <- rwmetrop(betabinexch,
                     list(var=laplace.fit$var, scale=2),
                     c(0, 0),
                     5000,
                     d[, c("Hits", "At.Bats")])


###################################################
### code chunk number 4: MultilevelModeling.Rnw:77-81
###################################################
mycontour(betabinexch, c(-1.5, -0.5, 2, 12), 
          d[, c("Hits", "At.Bats")],
          xlab="Logit ETA", ylab="Log K")
with(mcmc.fit, points(par))


###################################################
### code chunk number 5: MultilevelModeling.Rnw:88-99
###################################################
eta <- with(mcmc.fit, exp(par[, 1]) / (1 + exp(par[, 1])))
K <- exp(mcmc.fit$par[, 2])
p.estimate <- function(j, eta, K){
  yj <- d[j, "Hits"]
  nj <- d[j, "At.Bats"]
  p.sim <- rbeta(5000, yj + K * eta, nj - yj + K * (1 - eta))
  quantile(p.sim, c(0.05, 0.50, 0.95))  
}
E <- t(sapply(1:18, p.estimate, eta, K))
rownames(E) <- d[, "Name"]
round(E, 3)


###################################################
### code chunk number 6: MultilevelModeling.Rnw:105-115
###################################################
plot(d$Hits / 45, E[, 2], pch=19,
    ylim=c(.15, .40),
    xlab="Observed AVG", ylab="True Probability",
    main="90 Percent Probability Intervals")
for (j in 1:18)
  lines(d$Hits[j] / 45 * c(1, 1), E[j, c(1, 3)])
abline(a=0, b=1, col="blue")
abline(h=mean(d$Hits) / 45, col="red")
legend("topleft", legend=c("Individual", "Combined"),
       lty=1, col=c("blue", "red"))


