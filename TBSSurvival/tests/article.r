# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
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
####### Transform Both Sides Model: A parametric approach
#######

library("TBSSurvival")
## IF YOU DONT HAVE THE LIBRARY, PLEASE DOWNLOAD IT FROM
## http://code.google.com/p/tbssurvival/
## AND INSTALL IT USING:
## R CMD install TBSSurvival_version.tar.gz
## OR INSIDE R:
## install.packages("TBSSurvival_version.tar.gz",repos=NULL,type="source")
##

## To perform convergence analysis of BE set flag.convergence "TRUE".
## This takes some extra time to run however...
flag.convergence <- TRUE

initial.time <- proc.time()
# Seting initial seed to have the same results that we obtained.
set.seed(1234)

##########################################
## Section 2. Transform-both-side model ##
##########################################

i.fig <- 1 # counter figure number
## Figure 1(a)
## hazard functions with normal N(0,1) error
## plotted for lambda = 0.1, 1 and 2, 
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
t <- seq(0.01,10,0.01)
plot(t,htbs(t,lambda=0.1,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=1,col="gray60",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,htbs(t,lambda=1.0,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=2,col="gray40")
lines(t,htbs(t,lambda=2.0,xi=sqrt(1),beta=1,dist="norm"),type="l",lwd=2,lty=3,col="gray20")
title(xlab="time",ylab="h(t)",main=expression(epsilon ~ textstyle(paste("~",sep="")) ~ Normal(0,sigma^2==1)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(lambda == 0.1),
                expression(lambda == 1.0),
                expression(lambda == 2.0)),
              col=c("gray60","gray40","gray20"),lty=c(1,2,3),cex=1.1,lwd=2,bg="white")
dev.off()

## Figure 1(b)
## hazard functions with double-exp error (b=1)
## and lambda = 0.5, 1, and 2
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
t <- seq(0.01,10,0.01)
plot(t,htbs(t,lambda=0.5,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=1,col="gray60",
     axes="FALSE",ylim=c(0,1),xlab="",ylab="")
lines(t,htbs(t,lambda=1.0,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=2,col="gray40")
lines(t,htbs(t,lambda=2.0,xi=sqrt(1),beta=1,dist="doubexp"),type="l",lwd=2,lty=3,col="gray20")
title(xlab="time",ylab="h(t)",main=expression(epsilon ~ textstyle(paste("~",sep="")) ~ doubexp(b==1)),cex.lab=1.2)
axis(1,at=seq(0,10,1),lwd=2,lwd.ticks=2,pos=0)
axis(2,at=seq(0,1,0.2),lwd=2,lwd.ticks=2,pos=0)
legend(6,0.95,c(expression(lambda == 0.5),
                expression(lambda == 1.0),
                expression(lambda == 2.0)),
              col=c("gray60","gray40","gray20"),lty=c(1,2,3),cex=1.1,lwd=2,bg="white")
dev.off()


################################
##  Section 4 - Experiments   ##
################################


################################
# Section 4.1 - Real data sets #
################################

############################################################################################
# Data                                                                                     #
# Number of Cycles (in Thousands) of Fatigue Life for 67 of 72 Alloy T7987 Specimens that  #
# Failed Before 300 Thousand Cycles.                                                       #
# Meeker & Escobar, pp. 130-131 (1998)                                                     #
############################################################################################
data(alloyT7987)

## Summary statistic of the failure time in thousand cycles for the Alloy T7987 data set.
cat("Time AlloyT7987: Summary statistics\n")
print(summary(alloyT7987$time))

###########
## Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

## MLE Estimation for each possible error distribution
tbs.mle <- tbs.survreg.mle(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("all"))

## Building the Table of Example (MLE)
## This is only for ease of inclusion in the paper, as the results are all available
## in the corresponding R variables
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Some quantities of the TBS Model for the Alloy T7987 data set ",
                        "using MLE: AIC, BIC, parameter estimates and their standard ",
                        "deviations in parenthesis.}\n",
              "\\label{table_alloy-mle}\n",
              "\\centering\n",
              "\\begin{tabular}{cccccc}\n",
              "\\hline\n",
              "Error Distribution & AIC & BIC & $\\widehat{\\lambda}$ & $\\widehat{\\beta_0}$ & $\\widehat{\\xi}$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",sprintf("%.2f",tbs.mle$norm$AIC),       " & ",sprintf("%.2f",tbs.mle$norm$BIC),        " & ",
                                 sprintf("%.4f",tbs.mle$norm$lambda),    " & ",sprintf("%.4f",tbs.mle$norm$beta),     " & ",sprintf("%.4f",tbs.mle$norm$xi),     " \\\\ \n",
              "& & & {\\footnotesize(",sprintf("%.4f",tbs.mle$norm$lambda.se),")} & {\\footnotesize(",sprintf("%.4f",tbs.mle$norm$beta.se),")} & {\\footnotesize(",sprintf("%.5f",tbs.mle$norm$xi.se),")} \\\\ \n",
              "DoubExp       & ",sprintf("%.2f",tbs.mle$doubexp$AIC),       " & ",sprintf("%.2f",tbs.mle$doubexp$BIC),        " & ",
                                 sprintf("%.4f",tbs.mle$doubexp$lambda),    " & ",sprintf("%.4f",tbs.mle$doubexp$beta),     " & ",sprintf("%.4f",tbs.mle$doubexp$xi),     " \\\\ \n",
              "& & & {\\footnotesize(",sprintf("%.4f",tbs.mle$doubexp$lambda.se),")} & {\\footnotesize(",sprintf("%.4f",tbs.mle$doubexp$beta.se),")} & {\\footnotesize(",sprintf("%.5f",tbs.mle$doubexp$xi.se),")} \\\\ \n",
              "t-Student     & ",sprintf("%.2f",tbs.mle$t$AIC),       " & ",sprintf("%.2f",tbs.mle$t$BIC),        " & ",
                                 sprintf("%.4f",tbs.mle$t$lambda),    " & ",sprintf("%.4f",tbs.mle$t$beta),     " & ",sprintf("%.4f",tbs.mle$t$xi),     " \\\\ \n",
              "& & & {\\footnotesize(",sprintf("%.4f",tbs.mle$t$lambda.se),")} & {\\footnotesize(",sprintf("%.4f",tbs.mle$t$beta.se),")} & {\\footnotesize(",sprintf("%.5f",tbs.mle$t$xi.se),")} \\\\ \n",
              "Cauchy        & ",sprintf("%.2f",tbs.mle$cauchy$AIC),       " & ",sprintf("%.2f",tbs.mle$cauchy$BIC),        " & ",
                                 sprintf("%.4f",tbs.mle$cauchy$lambda),    " & ",sprintf("%.4f",tbs.mle$cauchy$beta),     " & ",sprintf("%.4f",tbs.mle$cauchy$xi),     " \\\\ \n",
              "& & & {\\footnotesize(",sprintf("%.4f",tbs.mle$cauchy$lambda.se),")} & {\\footnotesize(",sprintf("%.4f",tbs.mle$cauchy$beta.se),")} & {\\footnotesize(",sprintf("%.5f",tbs.mle$cauchy$xi.se),")} \\\\ \n",
              "Logistic      & ",sprintf("%.2f",tbs.mle$logistic$AIC),       " & ",sprintf("%.2f",tbs.mle$logistic$BIC),        " & ",
                                 sprintf("%.4f",tbs.mle$logistic$lambda),    " & ",sprintf("%.4f",tbs.mle$logistic$beta),     " & ",sprintf("%.4f",tbs.mle$logistic$xi),     " \\\\ \n",
              "& & & {\\footnotesize(",sprintf("%.4f",tbs.mle$logistic$lambda.se),")} & {\\footnotesize(",sprintf("%.4f",tbs.mle$logistic$beta.se),")} & {\\footnotesize(",sprintf("%.5f",tbs.mle$logistic$xi.se),")} \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_alloy-mle.tex")
rm(text)

## Estimating the Kaplan-Meier
## This is used to plot the graph with the comparison of the estimations
km <- survival::survfit(formula = survival::Surv(alloyT7987$time, alloyT7987$delta == 1) ~ 1)

i.fig <- i.fig+1
## Figure 2(a) - Reliability functions (MLE)
## Kaplan-Meier is also plotted
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(tbs.mle$norm,lwd=2,col="gray20",ylab="R(t)",
     xlab="t: number of cycles (in thousands)",
     main="Reliability function (MLE)",cex.lab=1.2)
km <- survival::survfit(formula = survival::Surv(alloyT7987$time, alloyT7987$delta == 1) ~ 1)
lines(km)
dev.off()

## Figure 2(b) - Hazard function (MLE)
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(tbs.mle$norm,plot.type="hazard",lwd=2,col="gray20",ylab="h(t)",
     xlab="t: number of cycles (in thousands)",
     main="Hazard function (MLE)",cex.lab=1.2)
dev.off()

## Summary of TBS model with normal error distribution (MLE).
summary(tbs.mle$norm)

## Quantile Estimation
cat("\n"); cat("Quatile estimates:\n")
cat("Median time: ",round(exp(tbs.mle$norm$beta),2),"\n",sep="")
cat("95% c.i.: (",round(exp(tbs.mle$norm$beta+qnorm(0.025,0,1)*tbs.mle$norm$beta.se),2),", ",
                  round(exp(tbs.mle$norm$beta+qnorm(0.975,0,1)*tbs.mle$norm$beta.se),2),")\n",sep="")

###########
## Part of the Bayesian estimation
cat("\n"); cat("Bayesian Estimation:\n")

## Bayesian Estimation
## We run for each of the five error functions
tbs.bayes.norm     <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("norm"),
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.t        <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("t"),
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.cauchy   <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("cauchy"),
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.doubexp  <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("doubexp"),
                                     burn=500000,jump=2000,size=1000,scale=0.07)
tbs.bayes.logistic <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("logistic"),
                                     burn=500000,jump=2000,size=1000,scale=0.07)

## Here it is the code for Building the Table that goes into the paper
## As before, all the info is already available in the R variables. The
## table is just for presentation purposes.
text <- paste("\\begin{center}\n",
              "\\begin{table}\n",
              "\\caption{Some quantities of the TBS Model for the Alloy T7987 data set using",
                        "BE: DIC, parameter estimates and their standard deviations in parenthesis.}\n",
              "\\label{table_alloy-bayes}\n",
              "\\centering\n",
              "\\begin{tabular}{ccccc}\n",
              "\\hline\n",
              "Error Distribution & DIC & $\\widehat{\\lambda}$ & $\\widehat{\\beta_0}$ & $\\widehat{\\xi}$ \\\\ \n",
              "\\hline \n",
              "\\hline \n",
              "Normal        & ",round(tbs.bayes.norm$DIC,2),       " & ",
                                 round(tbs.bayes.norm$lambda,4),    " & ",round(tbs.bayes.norm$beta,4),     " & ",round(tbs.bayes.norm$xi,4),     " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.norm$lambda.sd,4),")} & {\\footnotesize(",round(tbs.bayes.norm$beta.sd,4),")} & {\\footnotesize(",round(tbs.bayes.norm$xi.sd,5),")} \\\\ \n",
              "DoubExp       & ",round(tbs.bayes.doubexp$DIC,2),    " & ",
                                 round(tbs.bayes.doubexp$lambda,4), " & ",round(tbs.bayes.doubexp$beta,4),  " & ",round(tbs.bayes.doubexp$xi,4),  " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.doubexp$lambda.sd,4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$beta.sd,4),")} & {\\footnotesize(",round(tbs.bayes.doubexp$xi.sd,5),")} \\\\ \n",
              "t-Student     & ",round(tbs.bayes.t$DIC,2),          " & ",
                                 round(tbs.bayes.t$lambda,4),       " & ",round(tbs.bayes.t$beta,4),        " & ",round(tbs.bayes.t$xi,4),        " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.t$lambda.sd,4),")} & {\\footnotesize(",round(tbs.bayes.t$beta.sd,4),")} & {\\footnotesize(",round(tbs.bayes.t$xi.sd,5),")} \\\\ \n",
              "Cauchy        & ",round(tbs.bayes.cauchy$DIC,2),     " & ",
                                 round(tbs.bayes.cauchy$lambda,4),  " & ",round(tbs.bayes.cauchy$beta,4),   " & ",round(tbs.bayes.cauchy$xi,4),   " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.cauchy$lambda.sd,4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$beta.sd,4),")} & {\\footnotesize(",round(tbs.bayes.cauchy$xi.sd,5),")} \\\\ \n",
              "Logistic      & ",round(tbs.bayes.logistic$DIC,2),   " & ",
                                 round(tbs.bayes.logistic$lambda,4)," & ",round(tbs.bayes.logistic$beta,4), " & ",round(tbs.bayes.logistic$xi,4), " \\\\ \n",
              "& & {\\footnotesize(",round(tbs.bayes.logistic$lambda.sd,4),")} & {\\footnotesize(",round(tbs.bayes.logistic$beta.sd,4),")} & {\\footnotesize(",round(tbs.bayes.logistic$xi.sd,5),")} \\\\ \n",
              "\\hline\n",
              "\\end{tabular}\n",
              "\\end{table}\n",
              "\\end{center}\n",
              sep="")
cat(text, file="table_alloy-bayes.tex")
rm(text)

i.fig <- i.fig+1
## Figure 3(a) - Reliability functions / Boundary (Bayes)
## This code generates the plot of the Kaplan-Meier and TBS BE estimation with
## the logistic error distribution. It also plots the 95% HPD credible interval
name <- paste("tbs_fig_",i.fig,"a",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(tbs.bayes.logistic)
lines(km)
dev.off()

## Figure 3(b) - Hazard functions (Bayes)
## Again, the hazard function plus the 95% HPD interval for the TBS BE with logistic error
name <- paste("tbs_fig_",i.fig,"b",".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(tbs.bayes.logistic,plot.type="hazard")
dev.off()

#########################
## Quantile Estimation ##

## Quantile Estimation
cat("\n"); cat("Quatile estimates:\n")
cat("Median time: ",round(median(exp(tbs.bayes.logistic$post[,3])),2),"\n")
hpd <- HPDinterval(as.mcmc(exp(tbs.bayes.logistic$post[,3])),0.95)
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")


## Convergence analysis of the chain
if (flag.convergence) {
  # To perform convergence analysis set "TRUE" in the if above.
  # See the first lines of this code.
  
  i.conv <- 1
  ## Figure - ACF and Time series graphics
  ## Although this generates graphs regarding the convergence analysis,
  ## they do not appear in the body of the paper.
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
  par(mfrow=c(2,3))
  plot(ts(tbs.bayes.logistic$post[,1]),xlab="iteration",ylab=expression(lambda),main="",type="l")
  plot(ts(tbs.bayes.logistic$post[,2]),xlab="iteration",ylab=expression(beta[0]),main="",type="l")
  plot(ts(tbs.bayes.logistic$post[,3]),xlab="iteration",ylab=expression(xi),main="",type="l")
  acf(tbs.bayes.logistic$post[,1],main=expression(lambda),ci.col="gray40")
  acf(tbs.bayes.logistic$post[,2],main=expression(beta[0]),ci.col="gray40")
  acf(tbs.bayes.logistic$post[,3],main=expression(xi),ci.col="gray40")
  par(mfrow=c(1,1))
  dev.off()
  
  # Analysis with 4 chains
  chain.logistic1 <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=tbs.mle$logistic$beta,
                                    guess.lambda=tbs.mle$logistic$lambda,
                                    guess.xi=tbs.mle$logistic$xi,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic2 <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=6,
                                    guess.lambda=1,
                                    guess.xi=1,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic3 <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=4,
                                    guess.lambda=2,
                                    guess.xi=2,
                                    burn=1,jump=1,size=250000,scale=0.06)
  chain.logistic4 <- tbs.survreg.be(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist="logistic",
                                    guess.beta=5.5,
                                    guess.lambda=0.5,
                                    guess.xi=0.5,
                                    burn=1,jump=1,size=250000,scale=0.06)
  
  n <- length(chain.logistic1$post[,1])
  aux <- seq(1,n,500) 
  ## Figure: Ergodic Means - lambda
  ## Again, the graphs about convergence are generated, but do not appear in the paper.
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
  aux.1 <- cumsum(chain.logistic1$post[,1])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,1])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,1])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,1])/seq(1,n,1)
  plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(lambda),main="",
       ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()

  #Figure: Ergodic Means - xi
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
  aux.1 <- cumsum(chain.logistic1$post[,2])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,2])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,2])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,2])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(xi),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()

  #Figure: Ergodic Means - beta
  i.conv <- i.conv+1
  name <- paste("tbs_fig-conv_",i.conv,".eps",sep="")
  postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
  aux.1 <- cumsum(chain.logistic1$post[,3])/seq(1,n,1)
  aux.2 <- cumsum(chain.logistic2$post[,3])/seq(1,n,1)
  aux.3 <- cumsum(chain.logistic3$post[,3])/seq(1,n,1)
  aux.4 <- cumsum(chain.logistic4$post[,3])/seq(1,n,1)
   plot(aux,aux.1[aux],type="l",xlab="iteration",ylab=expression(beta[0]),main="",
        ylim=c(min(aux.1,aux.2,aux.3,aux.4),max(aux.1,aux.2,aux.3,aux.4)))
  lines(aux,aux.2[aux],type="l",col="gray30")
  lines(aux,aux.3[aux],type="l",col="gray20")
  lines(aux,aux.4[aux],type="l",col="gray10")
  dev.off()
  rm(aux,aux.1,aux.2,aux.3,aux.4)
  
  ## Gelman and Rubin Statistics
  ## It is presented in the paper to mention the convergence
  ## of the parameters estimated with BE
  R <- rep(0,3)
  for (j in 1:3) {
    W <- (var(chain.logistic1$post[,j])+var(chain.logistic2$post[,j])+var(chain.logistic3$post[,j])+var(chain.logistic4$post[,j]))/4
    mean.theta1 <- mean(chain.logistic1$post[,j])
    mean.theta2 <- mean(chain.logistic2$post[,j])
    mean.theta3 <- mean(chain.logistic3$post[,j])
    mean.theta4 <- mean(chain.logistic4$post[,j])
    mean.theta <- mean(c(chain.logistic1$post[,j],chain.logistic2$post[,j],chain.logistic3$post[,j],chain.logistic4$post[,j]))
    B <- (n/3)*((mean.theta1-mean.theta)^2+(mean.theta2-mean.theta)^2+(mean.theta3-mean.theta)^2+(mean.theta4-mean.theta)^2)
    R[j] <- (((1-1/n)*W+B/n)/W)
  }

  cat("\n"); cat("Gelman-Rubin Statistics:\n")
  cat("lambda: ",round(R[1],3),"\n")
  cat("    xi: ",round(R[2],3),"\n")
  cat("  beta: ",round(R[3],3),"\n")
}

############################################################################################
# Data                                                                                     #
# Colon Cancer - library(survival)                                                         #
############################################################################################
cat("\n"); cat("\n"); 
cat("++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("Data - Colon\n")
data(colon)

###########
# Part MLE
cat("\n"); cat("Maximum Likelihood Estimation:\n")

# Estimating TBS model with normal error distribution
fit.mle <- tbs.survreg.mle(survival::Surv(colon$time,colon$status) ~ colon$node4,dist=dist.error("norm"))

# Kaplan-Meier estimates
km <- survival::survfit(formula = survival::Surv(colon$time,colon$status) ~ colon$node4)

## Figure 4(a): estimated survival functions (MLE)
## This plots the Kaplan-Meier and the TBS curves for the colon data set,
## with the patients split into 2 groups based on the variable node4
i.fig <- i.fig+1
name <- paste("tbs_fig_",i.fig,".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(km,xlim=c(0,3000),conf.int=FALSE,lty=1,lwd=1,mark.time=FALSE,xlab="",ylab="")
title(ylab="S(t)",xlab="time",main="Survival function (MLE)",cex.lab=1.2)
lines(fit.mle,lty=c(1,2),lwd=c(3,3),col=c("gray50","gray50"))
dev.off()

## Quantile Estimation
## This results of the quantile estimation with MLE are not displayed in the paper,
## because they are very similar to the BE results, which we present.
cat("\n"); cat("Quatile estimates:\n")
cat("Median time | X=0: ",sprintf("%.2f",exp(fit.mle$beta[1])),"\n")
cat("95% c.i.: (",sprintf("%.2f",exp(fit.mle$beta[1]+qnorm(0.025,0,1)*fit.mle$beta.se[1])),",",
                  sprintf("%.2f",exp(fit.mle$beta[1]+qnorm(0.975,0,1)*fit.mle$beta.se[1])),")\n")
cat("\n")
cat("Median time | X=1: ",sprintf("%.2f",exp(sum(fit.mle$beta))),"\n")
cat("95% c.i.: (",sprintf("%.2f",exp(sum(fit.mle$beta)+qnorm(0.025,0,1)*sqrt(sum(fit.mle$beta.se^2)))),",",
                  sprintf("%.2f",exp(sum(fit.mle$beta)+qnorm(0.975,0,1)*sqrt(sum(fit.mle$beta.se^2)))),")\n")
cat("\n")
cat("Odds Median: ",sprintf("%.2f",exp(fit.mle$beta[2])),"\n")
cat("95% c.i.: (",sprintf("%.2f",exp(fit.mle$beta[2]+qnorm(0.025,0,1)*fit.mle$beta.se[2])),",",
                  sprintf("%.2f",exp(fit.mle$beta[2]+qnorm(0.975,0,1)*fit.mle$beta.se[2])),")\n")


###########
# Part BE
cat("\n"); cat("Bayesian Estimation:\n")

## Estimating TBS model with normal error distribution
## The same as before for the colon data set with node4 as discriminative covariate,
## but now using the Bayesian estimation
fit.be <- tbs.survreg.be(survival::Surv(colon$time,colon$status) ~ colon$node4,dist="norm",
                         burn=50000,jump=500,size=1000,scale=0.05)

## Figure 4(b): estimated survival functions (BE)
## It also shows the Kaplan-Meier and the 95% HPD credible intervals
i.fig <- i.fig+1
name <- paste("tbs_fig_",i.fig,".eps",sep="")
postscript(name,width=5,height=5,paper="special",colormodel="gray",horizontal=FALSE)
plot(km,xlim=c(0,3000),conf.int=FALSE,lty=1,lwd=1,mark.time=FALSE,xlab="",ylab="")
title(ylab="S(t)",xlab="time",main="Survival function (BE)",cex.lab=1.2)
lines(fit.be,lty=c(1,2),lwd=c(3,3),col=c("gray50","gray50"),
      lty.HPD=c(3,3),lwd.HPD=c(2,2),col.HPD=c("gray50","gray50"))
dev.off()

## Quantile Estimation
## This results are discussed in the final part of section 4.1 of the paper
cat("\n"); cat("Quatile estimates:\n")
median.0 <- mean(exp(fit.be$post[,3]))
hpd <- HPDinterval(as.mcmc(exp(fit.be$post[,3])),0.95)
cat("Median time | X=0: ",round(median.0,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
median.1 <- mean(exp(apply(fit.be$post[,3:4],1,sum)))
hpd <- HPDinterval(as.mcmc(exp(apply(fit.be$post[,3:4],1,sum))),0.95)
cat("Median time | X=1: ",round(median.1,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
O <- mean(exp(fit.be$post[,4]))
hpd <- HPDinterval(as.mcmc(exp(fit.be$post[,4])),0.95)
cat("Odds Median: ",round(O,2),"\n")
cat("95% HPD c.i.: (",round(hpd[1],2),",",round(hpd[2],2),")\n")
rm(median.0,median.1,O)


#####################
# End

## finally, just print the total time to run it.
## Without running the article_simulation.r (which many take a very long time),
## this file should run in about 2 to 3 hours in a modern computer (as of June-2012).
run.time <- (proc.time()[3]-initial.time[3])/60
print(run.time)
