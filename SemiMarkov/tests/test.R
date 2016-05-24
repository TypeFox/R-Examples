library(numDeriv)
library(MASS)
library(SemiMarkov)


## example with covariates for 3-state model
## Asthma control data
data(asthma)

## Defining vector of states and matrix of possible transitions
states.asthma <- c("1","2","3")
mtrans.asthma <- matrix(FALSE, ncol = 3, nrow = 3)
mtrans.asthma[1, 2:3] <- c("W","EW")
mtrans.asthma[2, c(1,3)] <- c("E","E")
mtrans.asthma[3, c(1,2)] <- c("W","W")

## checking the transitions in dataset
## censorship defined as "i->i", i=1,2,3
table.state(data = asthma, states = states.asthma, mtrans = mtrans.asthma)

## covariates
BMI <- as.data.frame(asthma$BMI)
SEX <- as.data.frame(asthma$Sex)

## semiMarkov model with Exponentiated Weibull waiting time distribution
## covariate BMI influencing transitions 1->2,1->3 and default initial values
sm.asthmaBMI2<-semiMarkov(data = asthma,cov =cbind( BMI,SEX), states = states.asthma, mtrans = mtrans.asthma, cov_tra=list(c("12","13"),c("21","23")))

sm.asthmaBMI2<-semiMarkov(data = asthma, states = states.asthma, mtrans = mtrans.asthma)

plot(hazard(sm.asthmaBMI2))

init <- param.init(data = asthma, cov = as.data.frame(asthma$BMI), states = states.asthma, mtrans = mtrans.asthma,cov_tra=list(c("12","32")),dist_init=c(rep(1.5,6),rep(1.8,4),rep(2,1)),proba_init=rep(0.05,6),coef_init=rep(2,2))




## Inputting states vector, model's transtions and censorship
states <- c("1","2","3")
trans <- matrix(FALSE, nrow = 3, ncol = 3)
trans[1, 2:3] <- c("W","W")
trans[2, c(1,3)] <- c("EW","W")
trans[3, c(1,2)] <- c("W","W")

## Default initial parameters of the model without covariates
## and Weibull as waiting time distribution
init_def <- param.init(data = asthma, states = states, 
              mtrans = trans)


## Initial model with values chosen
## Covariate "dage" influencing transitions " 1-> 2" and "3 -> 2"
init <- param.init(data = asthma, cov = as.data.frame(asthma$Sex),
    states = states, mtrans = trans, 
    cov_tra=list(c("12","32")),dist_init=c(rep(1.5,6),rep(1.8,6),rep(2,1)),
    proba_init=c(0.2,0.8,0.3,0.7,0.35,0.65),coef_init=rep(0.2,2))

##Usage of param.init for a 3-states death-illness model without dataset
## Defining vector of states
## 1 - healthy, 2 - ill, 3 - dead
states_th <- c("1","2","3")
mtrans_th <- matrix(FALSE, nrow = 3, ncol = 3)
mtrans_th[1,c(2,3)]<-c("E","W")
mtrans_th[2,c(1,3)]<-"EW"

init_th<-param.init(states=states_th, mtrans=mtrans_th, 
    proba_init=c(0.7,0.3,0.2,0.8))

