### R code from vignette source 'BayesFactors.Rnw'

###################################################
### code chunk number 1: BayesFactors.Rnw:34-46
###################################################
fire.counts <- c(75, 88, 84, 99, 79, 68, 86, 109, 73, 85, 101, 85,
                 75, 81, 64, 77, 83, 83, 88, 83, 78, 83, 78, 80,
                 82, 90, 74, 72, 69, 72, 76, 76, 104, 86, 92, 88)
hist(fire.counts, probability=TRUE, ylim=c(0, .08))
x <- 60:110
lines(x, dpois(x, lambda=mean(fire.counts)), col="red")
lines(x, dnorm(x, mean=mean(fire.counts), sd=12), col="blue")
lines(x, dnorm(x, mean=mean(fire.counts), sd=6), col="green")
legend("topright", legend=c("M1: Poisson(theta)",
                            "M2: N(theta, 12)",
                            "M3: N(theta, 6)"),
       col=c("red", "blue", "green"), lty=1)


###################################################
### code chunk number 2: BayesFactors.Rnw:67-71
###################################################
model.1 <- function(theta, y){
  sum(log(dpois(y, theta))) + 
    dgamma(theta, shape=280, rate=4)
}


###################################################
### code chunk number 3: BayesFactors.Rnw:74-77
###################################################
library(LearnBayes)
log.pred.1 <- laplace(model.1, 80, fire.counts)$int
log.pred.1


###################################################
### code chunk number 4: BayesFactors.Rnw:81-91
###################################################
model.2 <- function(theta, y){
  sum(log(dnorm(y, theta, 6))) + 
    dgamma(theta, shape=280, rate=4)
}
model.3 <- function(theta, y){
  sum(log(dnorm(y, theta, 12))) + 
    dgamma(theta, shape=280, rate=4)
}
log.pred.2 <- laplace(model.2, 80, fire.counts)$int
log.pred.3 <- laplace(model.3, 80, fire.counts)$int


###################################################
### code chunk number 5: BayesFactors.Rnw:95-96
###################################################
data.frame(Model=1:3, log.pred=c(log.pred.1, log.pred.2, log.pred.3))


###################################################
### code chunk number 6: BayesFactors.Rnw:99-100
###################################################
exp(log.pred.1 - log.pred.3)


