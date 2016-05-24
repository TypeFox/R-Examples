# IDPSurvival package for R (http://www.R-project.org)
# Copyright (C) 2014 Francesca Mangili, Alessio Benavoli, Marco Zaffalon, 
#                    Cassio de Campos.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######
####### Code to generate the results in the paper:
####### Imprecise Dirichlet Process for the estimate and comparison of 
####### survival functions with censored data
#######

## IF YOU DONT HAVE THE LIBRARY, PLEASE DOWNLOAD IT FROM
## ?????
## AND INSTALL IT USING:
## R CMD install IDPSurvival.tar.gz
## OR INSIDE R:
## install.packages("IDPSurvival.tar.gz",repos=NULL,type="source")
##

library(IDPsurvival)

################################################
##       Section 3 - sIMULATIONS              ##
################################################


exp.survsamp <- function(n, lambda) {
  deaths <- matrix(rexp(n, rate = lambda),n,1)
  censor <- matrix(rexp(n, rate = lambda),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}

## Code used for the results in

## Figure 1
## IDP and Kaplan-Meier survival curves for a sample of both survival and censoring 
## times drawn from the exponential distribution with hazard lambda=5. Plotted for 
## s = 0.25 and n=25 and 100.

n <- 100
lambda <- 5
s <- 0.25
alpha <- 0.05

samp20 <- exp.survsamp(n,lambda)
iout <- isurvfit(Surv(time,status)~1, samp20, conf.int= 1-alpha, s=s, display=TRUE)
out <- survfit(Surv(time,status)~1, samp20, conf.int= 1-alpha,conf.type="log-log")
lines(out, col='red',mark=3)
t<-seq(from = 0., to = 0.6, by = 0.01)
lines(t, exp(-lambda*t), col='green')
legend('bottomleft',c("IDP upper","IDP lower", "IDP cred. int.",
                      "KM","KM conf. int.",expression(exp(-lambda*t))),
       lty=c(1,1,2,1,2,1),lwd=c(2,1,1,1,1,1), 
       col=c('black','black','black','red','red','green'),
       pch=c('o','o','o','+','.','.'), text.width=0.04,inset=0,bty="n")


## Code used for the results in

## Figure 2(a)
## Coverage of the Kaplan-Meier estimator with plain and log-log confidence types 
## and of the IDP estimator with s = 0.25. Evaluated for n=25.

## Figure 2(b)
## Coverage of the Kaplan-Meier estimator with plain and log-log confidence types 
## and of the IDP estimator with s= 0.25. Evaluated for n=100.

## Table 1
## Coverage of the Kaplan-Meier estimator with plain and log-log confidence types 
## and of the IDP estimator with s= 0.01, 0.25, 0.5. 
## Evaluated for n=25 and 100. 

runs <- 1000
n.eval <- 10
step <- 0.05
n <- 100
lambda <- 5
s <- 0.5
t.eval <-matrix(seq(from = 0.025, to = step*n.eval, by = step),n.eval,1)
t.mat <- repmat(t.eval,1,n)
ref <- t(exp(-lambda*t.eval))
error <- matrix(numeric(0),n.eval,0)
llerror <- matrix(numeric(0),n.eval,0)
ierror <- matrix(numeric(0),n.eval,0)
sc<-0
for ( k in 1:runs) {
  if (k%%100 ==1) cat(k,"\n")
  samp <- exp.survsamp(n,lambda)
  
  # Estimate coverage for the Kaplan-Meier "plain" confidence interval
  out <- survfit(Surv(time,status)~1, samp, conf.int= 1-alpha,conf.type="log")
  time <- matrix(out$time,1,n)
  mat.time <- repmat(time,n.eval,1)
  time.ref <- apply((mat.time<t.mat),1,sum)+1
  # Set estimates before the first observed death to NA (not considered in the 
  # coverage estimation)
  n.death <- cumsum(out$n.event)
  out$lower[which(n.death==0)] <- NA
  out$upper[which(n.death==0)] <- NA
  # Add estimate in 0
  lower <- c(NA,out$lower)
  upper <-c(NA,out$upper)
  # 
  kerr <- matrix((lower[time.ref]<ref)*(upper[time.ref]>ref),n.eval,1)
  error <- cbind(error, kerr)
  
  # Estimate coverage for the the Kaplan-Meier "log-log" confidence interval
  out <- survfit(Surv(time,status)~1, samp, conf.int= 1-alpha,conf.type="log-log")
  time<-matrix(out$time,1,n)
  mat.time <- repmat(time,n.eval,1)
  time.ref <- apply((mat.time<t.mat),1,sum)+1
  n.death <- cumsum(out$n.event)
  # Set estimates before the first observed death to NA (not considered in the 
  # coverage estimation)
  out$lower[which(n.death==0)] <- NA
  out$upper[which(n.death==0)] <- NA
  # Add estimate in 0
  lower <- c(NA,out$lower)
  upper <-c(NA,out$upper)
  kerr <- matrix((lower[time.ref]<ref)*(upper[time.ref]>ref),n.eval,1)
  llerror <- cbind(llerror, kerr)
  
  # Estimate coverage for the IDP credible interval
  iout <- isurvfit(Surv(time,status)~1, samp, conf.int= 1-alpha, s=s, display=FALSE,nsamples=1000)
  ilower <- c(iout$lower0,iout$lower)
  iupper <-c(1, iout$upper)
  ikerr <- matrix((ilower[time.ref]<ref)*(iupper[time.ref]>ref),n.eval,1)
  ierror <- cbind(ierror, ikerr)
}

cov <- matrix(numeric(0),3,0) 
sumall <- rep(0,3)
nall <-0
cat("\n")
for (i in 1:n.eval) {
  # Remove NA (i.e., samples for which t.eval(i) is too small or too large)
  if (sum(is.na(error[i,]))>0) {
    i.err <- cbind(error[i,-which(is.na(error[i,]))],
                   llerror[i,-which(is.na(error[i,]))],
                   ierror[i,-which(is.na(error[i,]))]
                  )
  } else {
    i.err <- cbind(error[i,], llerror[i,], ierror[i,])
  }
  cov <- cbind(cov,apply(i.err,2,mean))
  sumall <- sumall+apply(i.err,2,sum)
  nall <- nall+dim(i.err)[1]
  cat(dim(i.err)[1],"\n")
}

plot(t.eval, cov[3,], xlab="Time", ylab="Coverage", type = "p",pch='x',ylim=c(0.5,1))
points(t.eval, cov[1,], type = "p")
points(t.eval, cov[2,], type = "p",pch='*')
lines(c(0,0.5),c(1-alpha,1-alpha))
legend('bottomleft',c("IDP","Kaplan-Meier (plain)","Kaplan Meier (log-log)"),pch=c('x','o','*'), border = NA)

cat("\nCoverage of the Kaplan-Meier confidence interval: ")
cat(sumall[1]/nall)
cat("\nCoverage of the Kaplan-Meier log-log confidence interval: ")
cat(sumall[2]/nall)
cat("\nCoverage of the IDP credible interval: ")
cat(sumall[3]/nall,"\n")


#################################################################################
## Code used for the results in 

## Example 1
## type I error for the log-rank, Peto-Peto and IDP (s=0.25) tests for two samples
## with size n1 = 20 and n2 = 100 drawn from and exponential distribution with
## hazard lambda = 1. Censoring times from uniform U[0,2]

## Table 2
## Power (or type I error) of the log-rank, Peto-Peto and IDP tests (with s = 0.25) 
## evaluated for the scenarios SH (same hazard) and PH (proportional hazard). 

maketest <- function(dataset,groups,group_label,s=0.25,display=TRUE,
                     alpha=0.05,nsamples=10000,smax=FALSE){
  colnames(dataset)[which(colnames(dataset)==group_label)]<-"lab"
  fdata <- Surv(time, status) ~ lab
  testLR <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 0)
  testPP <- survdiff(fdata, dataset, lab==groups[1]|lab==groups[2], rho = 1)
  Xlab <- paste('lab',groups[1],sep='=')
  ix <- which(as.data.frame(testLR$n)[,1]==Xlab)
  zLR <- sign(testLR$obs[ix]-testLR$exp[ix])*sqrt(testLR$chisq)
  pLR <- 1-pnorm(zLR)
  pLR2 <- 2*pnorm(-sqrt(testLR$chisq))
  hLR <- 1*(pLR<alpha)
  zPP <- sign(testPP$obs[ix]-testPP$exp[ix])*sqrt(testPP$chisq)
  pPP <- 1-pnorm(zPP)
  pPP2 <- 2*pnorm(-sqrt(testPP$chisq))
  hPP <- 1*(pPP<alpha)
  out <- isurvdiff(fdata,dataset,groups=groups, alternative = 'greater',
                   display=display, s=s, nsamples=nsamples)
  hIDP <- out$h
  
  if (display==TRUE) {
    out2 <- isurvfit(fdata,dataset,subset=(lab==groups[1]|lab==groups[2]),s=s)
    cat("Number of data in group 1/2: ",  out2$strata[1],"/ ",out2$strata[2],"\n")
    cat("Lower/upper probability:     ", out$prob[1],"/ ",out$prob[2], "\n")
    cat("Log-rank test (p-value):     ", pLR, "\n")
    cat("Peto & Peto test (p-value):  ", pPP, "\n")
  }
  
  if (smax==TRUE) {
    s.out <- isurvdiff.smax(fdata,dataset,groups=groups, alternative = 'greater',
                      nsamples=nsamples)
    cat("Min s giving indecision:     ", s.out, "\n")
  }
  
  temp <- c(hLR,hPP,hIDP,s.out) 
  temp
}


expunif.survsamp <- function(n, lambda, UB) {
  deaths <- matrix(rexp(n, rate = lambda),n,1)
  censor <- matrix(runif(n,UB[1],UB[2]),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}


runs <- 1
s <- 0.25
alpha <- 0.05
## Example 1
n1 <- 100
n2 <- 100
lambda1 <- 1
lambda2 <- 1
# ## Table 2 SH
# n1 <- 30
# n2 <- 30
# lambda1 <- 1
# lambda2 <- 1
# ## Table 2 PH
# n1 <- 30
# n2 <- 30
# lambda1 <- 2
# lambda2 <- 1

UB1 <- c(0,2)
UB2 <- c(0,2)
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
    if (k%%50 ==1) cat(k,"\n")
    Z1 <- expunif.survsamp(n1, lambda1, UB1)
    Z1$lab <- rep(0,n1)
    Z2 <- expunif.survsamp(n2, lambda2, UB2)
    Z2$lab <- rep(1,n2)
    dataset <- rbind(Z1,Z2)
    Htest[k,] <- maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
}

errI <- apply(1*(Htest==1),2,mean)
indecision <- mean(1*(Htest[,3]==2))

cat("Type I error\n")
cat("Log-rank test : ", errI[1], "\n")
cat("Peto-Peto test: ", errI[2], "\n")
cat("IDP test      : ", errI[3], " (", indecision, ")","\n")


####################################################################################
## Artificial scenarios used for the results of:

## Table 2
## Power (or type I error) of the log-rank, Peto-Peto and IDP tests (with s = 0.25) 
## evaluated for the scenarios EHD (early hazard difference) and LHD (late hazard 
## difference). 

## Figure 3(a)
## survival function for the two curves considered in the scenario A.

## Figure 3(b) 
## survival function for the two curves considered in the scenario B.

## Table 3
## Type I error of of the log-rank, Peto-Peto (alpha=0.01 and 0.05) and IDP 
## (s=0.25 and alpha=0.05) tests for scenarios A and B.

step.hazard.survsamp <- function(n, lvec, tvec, UB) {
  tvec <- c(0,tvec)
  dt <- diff(tvec)
  dhaz <- lvec[-length(lvec)]*dt
  chaz <- c(0,cumsum(dhaz))
  chaz.samp <- -log(matrix(runif(n,0,1),n,1))
  tbin <- as.numeric(cut(chaz.samp,chaz))
  tbin[is.na(tbin)]<-length(lvec)
  deaths <- (chaz.samp-chaz[tbin])/lvec[tbin]+tvec[tbin]
  censor <- matrix(runif(n,UB[1],UB[2]),n,1)
  status <- matrix((deaths<censor)*1,n,1)
  time <- deaths*status+censor*(1-status)
  temp <- data.frame(time,status)
}

runs <- 1000
s <- 0.25
alpha <- 0.025
## Early Hazard Difference (EHD) scenario
n1 <- 30
n2 <- 30
tvec <- c(0.4, 0.6)
lvec1 <- c(3, 0.75, 1)
lvec2 <- c(0.75, 3, 1)
# ## Late Hazard Difference (LHD) scenario
# n1 <- 30
# n2 <- 30
# tvec <- c(0.4)
# lvec1 <- c(2, 4)
# lvec2 <- c(2, 0.4)
# ## Case A scenario
# n1 <- 50
# n2 <- 25
# tvec <- c(0.4)
# lvec1 <- c(0.75, 5)
# lvec2 <- c(2, 0.25)
# ## Case B scenario
# n1 <- 25
# n2 <- 50
# tvec <- c(0.4)
# lvec1 <- c(2, 0.2)
# lvec2 <- c(0.75, 5)

UB <- c(0,2)
Htest <- matrix(numeric(0),runs,3)
for (k in 1:runs) {
  if (k%%50 ==1) cat(k,"\n")
  Z1 <- step.hazard.survsamp(n1, lvec1, tvec, UB) 
  Z1$lab <- rep(0,n1)
  Z2 <- step.hazard.survsamp(n2, lvec2, tvec, UB) 
  Z2$lab <- rep(1,n2)
  dataset <- rbind(Z1,Z2)
  Htest[k,] <- maketest(dataset,groups=c(0,1),group_label="lab",s=s,display=FALSE,alpha=alpha)
}

errI <- apply(1*(Htest==1),2,mean)
indecision <- mean(1*(Htest[,3]==2))

cat("Type I error\n")
cat("Log-rank test : ", errI[1], "\n")
cat("Peto-Peto test: ", errI[2], "\n")
cat("IDP test      : ", errI[3], " (", indecision, ")","\n")


###############################################################
##          Section 4 Australian AIDS survival data          ##
###############################################################


## Table 4
## p-values (Log-rank and Peto-Peto tests), posterior probabilities an maximum s 
## giving a determinate decision (IDP test ## with s=0.25) for the AIDS dataset. 

## Figures from 4
## Survival curves (s = 0.25) for the comparisons in Table 4.


s<-0.25

data(Aids2,package='MASS')
dataset <- Aids2
dataset["time"]<-dataset[4]-dataset[3]
dataset["agg"] <- 1*(dataset["T.categ"]=="hs")+1*(dataset["T.categ"]=="het")+
                  2*(dataset["T.categ"]=="hsid")+2*(dataset["T.categ"]=="id")  
dataset[5]<-as.numeric(unlist(dataset[5]))

cat("\nMale < Female\n")
t1 <- maketest(dataset,c("M","F"),"sex",s=s,smax=TRUE)
legend('bottomleft',c("Female (X\')","Male (X)"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")

cat("\nNSW < VIC \n")
t2 <- maketest(dataset,c("NSW","VIC"),"state",s=s,smax=TRUE)
legend('bottomleft',c("New South Wales (X)","Victoria (X\')"),
       lty=c(1,1),lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")

cat("\nQLD < NSW \n")
t3 <- maketest(dataset,c("QLD","NSW"),"state",s=s,smax=TRUE)
legend('bottomleft',c("New South Wales (X\')","Queensland (X)"),
       lty=c(1,1),lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")

cat("\nno drug user < drug user \n")
t4 <- maketest(dataset,c(1,2),"agg",s=s,smax=TRUE)
legend('bottomleft',c("No drug user (X)","Drug user (X\')"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")

cat("\nBlood < Haemophilia \n")
t5 <- maketest(dataset,c("blood","haem"),"T.categ",s=s,smax=TRUE)
legend('bottomleft',c("Hemophilia (X\')","Blood (X)"),lty=c(1,1),
       lwd=c(2,2),col=c('black','red'),pch=c(1,3),bty="n")


###################################################################################
# End


