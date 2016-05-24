## download and install the HiddenMarkov Package:
## http://homepages.paradise.net.nz/david.harte/SSLib/

## Once its done load STAR and HiddenMarkov
require(STAR)
require(HiddenMarkov)
originalpar <- par(ask = dev.interactive(orNone = TRUE))
oldpar <- par()

## Start with the last simulated data set in
## renewalTestPlot example
isi.nb <- 500
myTransition <- matrix(c(0.9,0.1,0.1,0.9),2,2,byrow=TRUE)
states2 <- numeric(isi.nb+1) + 1
for (i in 1:isi.nb) states2[i+1] <- rbinom(1,1,prob=1-myTransition[states2[i],])+1
myLnormPara2 <- matrix(c(log(0.01),0.25,log(0.05),0.25),2,2,byrow=TRUE)
train2 <-
cumsum(rlnorm(isi.nb+1,myLnormPara2[states2,1],myLnormPara2[states2,2]))
## make the renewal test
renewalTestPlot(train2);layout(matrix(1,1,1))

## built the histogram of log(ISI)
hist(log(diff(train2)),breaks=20)
## The data are nice and the 2 peaks can be nicely
## separated based on ISI < -3.75
abline(v=-3.75,lty=2)
## We therefore try a HHM with 2 components.
## We get an intial guess for the 2 components from
## the roughly split data set
comp1 <- c(mean(log(diff(train2)[log(diff(train2)) < -3.75])),
           sd(log(diff(train2)[log(diff(train2)) < -3.75])))
comp2 <- c(mean(log(diff(train2)[log(diff(train2)) >= -3.75])),
           sd(log(diff(train2)[log(diff(train2)) >= -3.75])))
## We then get a rough idea of how many ISIs are in each
## component using the "empirical cumulative distribution function
plot(ecdf(log(diff(train2))))
## We have an approximate 50/50 repartition
abline(h=0.51,lty=2,col=2)
## We get an initial guess for the transition matrix
## corresponding to this naive repartition
Pi <- matrix(0.5,nrow=2,ncol=2)
## after CHECKING the documentation of the Baum.Welch
## function we fit the data
fit1 <- Baum.Welch(x=log(diff(train2)),
                   Pi=Pi,
                   delta=c(1,0),
                   distn="norm",
                   pm=list(mean=c(comp1[1],comp2[1]),sd=c(comp1[2],comp2[2]))
                   )
## We can look at the estimated transition matrix
fit1$Pi
## and compare it with the true one
abs(myTransition-fit1$Pi)/myTransition
## We can look at the estiamted distribution para
fit1$pm
## and compare with "truth"
abs(unlist(fit1$pm)-as.numeric(myLnormPara2))/abs(as.numeric(myLnormPara2))

## example using neuron 2 of data set CAL2S
data(CAL2S)
CAL2S <- lapply(CAL2S,as.spikeTrain)
CAL2S[["neuron 2"]]
renewalTestPlot(CAL2S[["neuron 2"]],lag.max=50);layout(matrix(1,1,1))
summary(CAL2S[["neuron 2"]])
hist(log(diff(CAL2S[["neuron 2"]])),breaks=20)
abline(v=-2.5,lty=2)
compA <- c(mean(log(diff(CAL2S[["neuron 2"]])[log(diff(CAL2S[["neuron 2"]])) < -2.5])),
           sd(log(diff(CAL2S[["neuron 2"]])[log(diff(CAL2S[["neuron 2"]])) < -2.5])))
compB <- c(mean(log(diff(CAL2S[["neuron 2"]])[log(diff(CAL2S[["neuron 2"]])) >= -2.5])),
           sd(log(diff(CAL2S[["neuron 2"]])[log(diff(CAL2S[["neuron 2"]])) >= -2.5])))
plot(ecdf(log(diff(CAL2S[["neuron 2"]]))))
abline(v=-2.5,lty=2)
abline(h=0.7,lty=2)
Pi <- matrix(rep(c(0.7,0.3),2),nrow=2,byrow=TRUE)
fitA <- Baum.Welch(x=log(diff(CAL2S[["neuron 2"]])),
                   Pi=Pi,
                   delta=c(1,0),
                   distn="norm",
                   pm=list(mean=c(compA[1],compB[1]),sd=c(compA[2],compB[2]))
                   )
## We don't know the truth here but we can make some
## rough checks to start with
## We redraw the ISI histogram on the probability scale
hist(log(diff(CAL2S[["neuron 2"]])),breaks=20,prob=TRUE)
## we add the component with the short ISI
xx <- seq(-6,0,0.01)
lines(xx,
      sum(fitA$u[,1]>0.5)/dim(fitA$u)[1]*dnorm(xx,fitA$pm$mean[1],fitA$pm$sd[1]),
      col=2
      )
## we add the component with the long ISI
lines(xx,
      sum(fitA$u[,1]<=0.5)/dim(fitA$u)[1]*dnorm(xx,fitA$pm$mean[2],fitA$pm$sd[2]),
      col=4
      )
## we add the sum of the 2 components
lines(xx,
      sum(fitA$u[,1]>0.5)/dim(fitA$u)[1]*dnorm(xx,fitA$pm$mean[1],fitA$pm$sd[1])+
      sum(fitA$u[,1]<=0.5)/dim(fitA$u)[1]*dnorm(xx,fitA$pm$mean[2],fitA$pm$sd[2]),
      col=3
      )
## a better check is obtained by rescaling the ISI
nISI <- length(CAL2S[["neuron 2"]])-1
isiN2l <- log(diff(CAL2S[["neuron 2"]]))
rISI <- sapply(1:nISI,
               function(idx) {
                 if (fitA$u[idx,1]>0.5) {
                   qexp(pnorm(isiN2l[idx],fitA$pm$mean[1],fitA$pm$sd[1]))
                 } else {
                   qexp(pnorm(isiN2l[idx],fitA$pm$mean[2],fitA$pm$sd[2]))
                 }
               }
               )
rescaledTrain <- as.spikeTrain(cumsum(rISI))
## Now check it out
rescaledTrain
renewalTestPlot(rescaledTrain,lag.max=50);layout(matrix(1,1,1))
plot(varianceTime(rescaledTrain),style="Ogata")
## Rather good

## If not sure about the number of components, fit several models
## (ie one with 2, one with 3, etc). Each model has 2*k+k*(k-1)
## parameters. Use the AIC -2*logLikelihood + 2*number of parameters
## and keep the model with the smallest AIC. The component $LL of the list
## returned by Baum.Welch contains the loglikelihood.
