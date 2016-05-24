###################################################
### Chap15Start
###################################################
library(mistat)


###################################################
### BayesBinomBeta
###################################################
library(LearnBayes)

Probs <- c(0.5, 0.5)

BetaPar1<- c(1, 1)

BetaPar2 <- c(15, 2)

Betapar <- rbind(BetaPar1, BetaPar2)

Data<- c(10, 2)

binomial.beta.mix(probs=Probs,
                  betapar=Betapar,
                  data=Data)


###################################################
### BayesPoissonGamma
###################################################
Probs <- c(0.5, 0.5)

GammaPar1 <- c(1, 1)  # Gamma parameters are expressed as 
# shape and rate
# scale is 1/rate
GammaPar2 <- c(15, 2) 

Gammapar <- rbind(GammaPar1, GammaPar2)

Data<- list(
  y=c(5),
  t=c(1))

poisson.gamma.mix(probs=Probs, 
                  gammapar=Gammapar, 
                  data=Data)


###################################################
### Chap15End
###################################################
rm(Betapar, Gammapar, BetaPar1, BetaPar2, 
   GammaPar1, GammaPar2, Data, Probs)
detach(package:LearnBayes)
detach(package:mistat)
