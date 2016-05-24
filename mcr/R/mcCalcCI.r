###############################################################################
##
## mcCalcCI.R
##
## Functions for computing confidence intervals.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Analytical Confidence Interval
#'
#' Calculate wald confidence intervals for intercept and slope
#' given point estimates and standard errors.
#'
#' @param b0 point estimate of intercept.
#' @param b1 point estimate of slope.
#' @param se.b0 standard error of intercept.
#' @param se.b1 standard error of slope.
#' @param n number of observations.
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @return 2x4 matrix of estimates and confidence intervals for intercept and slope.
mc.analytical.ci <- function(b0,b1,se.b0,se.b1,n,alpha) 
{
	rmat <- matrix(nrow=2,ncol=4)
	rownames(rmat) <- c("Intercept","Slope")
	colnames(rmat) <- c("EST","SE","LCI","UCI")
	rmat[,1] <- c(b0,b1)
	rmat[,2] <- c(se.b0,se.b1)
	rmat["Intercept","LCI"] <- b0-qt(1-alpha/2,n-2)*se.b0
	rmat["Intercept","UCI"] <- b0+qt(1-alpha/2,n-2)*se.b0
	rmat["Slope","LCI"] <- b1-qt(1-alpha/2,n-2)*se.b1
	rmat["Slope","UCI"] <- b1+qt(1-alpha/2,n-2)*se.b1
	return(rmat)
}

#' Jackknife Confidence Interval
#'
#' Calculate Jackknife confidence intervals for intercept, slope or bias
#' given of vector of jackknife point estimates and global point estimate.
#'
#' @param Xjack vector of point estimates for jackknife samples. 
#'        The i-th element contains point estimate for the dataset without the i-th observation. 
#' @param xhat point estimate for the complete data set (scalar).
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @return a list with elements
#'  \item{est}{point estimate for the complete data set (scalar).}
#'  \item{se}{standard deviation of point estimate calculated with Jackknife Method.}
#'  \item{CI}{confidence interval for point estimate.}
#' @references  Linnet, K. (1993)
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              \emph{CLIN. CHEM.} \bold{39/3}, 424--432.
mc.calcLinnetCI<-function(Xjack,xhat,alpha=0.05) 
{
    npoints<-length(Xjack)
    delta.X <- npoints*xhat-(npoints-1)*Xjack
    se <- sd(delta.X)/sqrt(npoints)
    CI <- c(xhat-qt(1-alpha/2,npoints-2)*se , xhat+qt(1-alpha/2,npoints-2)*se)
    return(list(est=xhat,se=se,CI=CI))
}

#' Student Method for Calculation of Resampling Confidence Intervals
#'
#' Calculate bootstrap confidence intervals for intercept, slope or bias
#' given a vector of bootstrap point estimates.
#'
#' @param Xboot vector of point estimates for each bootstrap sample. 
#'        The i-th element contains the point estimate of the i-th bootstrap sample. 
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @param xhat global point estimate for which the confidence interval shall be computed.
#' @param npoints number of points used for the regression analysis.
#' @return a list with elements
#'  \item{est}{the point estimate xhat}
#'  \item{se}{standard deviation computed from bootstrap point estimates Xboot} 
#'  \item{CI}{Confidence interval for point estimate xhat, calculated as \eqn{xhat +/- qt(1-alpha,n-2)*sd}.}
#' @references  Carpenter, J., Bithell, J. (2000)
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              \emph{Stat Med}, \bold{19 (9)}, 1141--1164.     
mc.calc.Student<-function(Xboot, xhat, alpha, npoints)
{
    est <-xhat
    se <- sd(Xboot)
    CI <- c(est-qt(1-alpha/2,npoints-2)*se, est+qt(1-alpha/2, npoints-2)*se)
    return(list(est=xhat, se=se, CI=CI))
}

#' Bootstrap-t Method for Calculation of Resampling Confidence Intervals
#'
#' Calculate resampling confidence intervals for intercept, slope or bias
#' with t-Boot method given a vector of bootstrap point estimates and a vector of bootstrap 
#' standard deviations.
#'
#' @param Xboot vector of point estimates for bootstrap sample. 
#'        The i-th element contains the point estimate for the i-th bootstrap sample. 
#' @param Sboot vector of standard deviations for each bootstrap sample. 
#'        It schould be estimated with any analytical method or nonparametric with nested bootstrap.  
#' @param xhat point estimate for the complete data set (scalar).
#' @param shat estimate of standard deviation for the complete data set (scalar).
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @return a list with elements
#'  \item{est}{point estimate for the complete data set (xhat).}
#'  \item{se}{estimate of standard deviation for the complete data set (shat).} 
#'  \item{CI}{confidence interval for the point estimate.}
#' @references  Carpenter, J., Bithell, J. (2000)
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              \emph{Stat Med}, \bold{19 (9)}, 1141--1164.    
mc.calc.tboot <- function(Xboot, Sboot, xhat, shat, alpha)
{
    tstar <- (Xboot-xhat)/Sboot
    if(is.na(shat)) 
        shat <- sd(Xboot)
    CI <- c(xhat-mc.calc.quant(X=tstar,alpha=1-alpha/2)*shat,xhat-mc.calc.quant(X=tstar,alpha=alpha/2)*shat)
    return(list(est=xhat,se=shat,CI=CI))
}

#' Quantile Calculation for BCa
#' 
#' We are using the R default (SAS (type=3) seems bugged) quantile calculation instead of the quantile function
#' described in Effron&Tibshirani.
#'
#' @param X numeric vector.
#' @param alpha probabilty
#' @return alpha-quantile of vector X.
mc.calc.quant <- function(X, alpha) 
{
    #B<-length(X)
    #OX<-X[order(X)]   # 0 - der erster Index?
    #Q<-alpha*(B+1)
    #if (round(Q)==Q){
    #        tetta.Q<-OX[Q]
    #    } else {                              # Problem: wenn a=0
    #        a<-round(Q); b<-round(Q)+1
    #        tetta.a<-OX[a]
    #        tetta.b<-OX[b]
    #        tetta.Q <- tetta.a+
    #           ((qnorm(Q/(B+1))-qnorm(a/(B+1)))/(qnorm(b/(B+1))
    #           -qnorm(a/(B+1))))*(tetta.b-tetta.a)
    #        }
    #
    #return(tetta.Q)
    return(quantile(X, alpha))
}

#' Bias Corrected and Accelerated Resampling Confidence Interval
#'
#' Calculate resampling BCa confidence intervals for intercept, slope or bias
#' given a vector of bootstrap and jackknife point estimates.
#'
#' @param Xboot vector of point estimates for bootstrap samples. 
#'        The i-th element contains point estimate of the i-th bootstrap sample.
#' @param Xjack vector of point estimates for jackknife samples. 
#'        The i-th element contains point estimate of the dataset without i-th observation. 
#' @param xhat point estimate for the complete data set (scalar).
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @return a list with elements
#'  \item{est}{point estimate for the complete data set (xhat).}
#'  \item{CI}{confidence interval for point estimate.}
#' @references  Carpenter, J., Bithell, J. (2000)
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              \emph{Stat Med}, \bold{19 (9)}, 1141--1164.       
mc.calc.bca <- function(Xboot,Xjack,xhat,alpha)
{
    ## Assertions
    stopifnot(!any(is.na(Xboot)))
    ## Number of bootstrap samples
    B <- length(Xboot)
    ## Calculate fraction of bootstrap samples with Xboot<xhat
    ## taking into account weird special cases
    p <- sum(Xboot<xhat)
    if(p==0) 
        p <- 0.5 # 0 would cause crash due to infinite qnorm
    if(p==B) 
        p <- B-0.5 # B would cause crash due to infinite qnorm
    if((length(unique(Xboot))==1) && (unique(Xboot)==xhat)) 
        p <- B/2 # everything is the same, uninformative prob estimate
    ## Compute b
    b <- qnorm(p/B)
    ## Compute a
    Xtilde <- rep(mean(Xjack), length(Xjack))
    if(sum((Xtilde-Xjack)^2)>0)
        a <- sum((Xtilde-Xjack)^3)/(6*(sum((Xtilde-Xjack)^2))^(3/2))
    else 
        a <- 0
    ## Compute quantiles
    zL <- qnorm(alpha/2)
    zU <- qnorm(1-alpha/2)
    IU <- pnorm(b-(zL-b)/(1+a*(zL-b)))
    IL <- pnorm(b-(zU-b)/(1+a*(zU-b)))    
    QL <- mc.calc.quant(Xboot,IL)
    QU <- mc.calc.quant(Xboot,IU)
    return(list(est=xhat,CI=c(QL,QU)))
}

#' Quantile Method for Calculation of Resampling Confidence Intervals
#'
#' Calculate bootstrap confidence intervals for intercept, slope or bias
#' given the vector of bootstrap point estimates.
#'
#' @param Xboot vector of point estimates for bootstrap samples. 
#'        The i-th element contains point estimate of the i-th bootstrap sample. 
#' @param alpha numeric value specifying the 100(1-alpha)\% confidence level for the confidence interval (Default is 0.05).
#' @return a list with elements
#'  \item{est}{median of bootstrap point estimates Xboot.} 
#'  \item{CI}{confidence interval for point estimate 'est', calculated as quantiles.}
#' @references  B. Efron and RJ. Tibshirani (1994)
#'              An Introduction to the Bootstrap.
#'              \emph{Chapman & Hall}.     
mc.calc.quantile <- function(Xboot, alpha) 
{
    stopifnot(!any(is.na(Xboot)))
    est <- quantile(Xboot, probs=0.5)
    CI <- quantile(Xboot, probs=c(alpha/2, 1-alpha/2))
    return(list(est=est,CI=CI))
}

#' Returns Results of Calculations in Matrix Form
#'
#' @param b0 point estimate for intercept.
#' @param b1 point estimate for slope.
#' @param se.b0 standard error of intercept estimate.
#' @param se.b1 standard error of slope estimate.
#' @param CI.b0 numeric vector of length 2 - confidence interval for intercept.
#' @param CI.b1 numeric vector of length 2 - confidence interval for slope.
#' @return 2x4 matrix of estimates and confidence intervals for intercept and slope.
mc.make.CIframe <- function(b0, b1, se.b0, se.b1, CI.b0, CI.b1) 
{
    bci <- matrix(NA, nrow=2, ncol=4)
    colnames(bci) <- c("EST","SE","LCI","UCI")
    rownames(bci) <- c("Intercept","Slope")
    bci[,"EST"] <- c(b0,b1)
    bci[,"SE"] <- c(se.b0, se.b1)
    bci["Intercept", c("LCI","UCI")] <- CI.b0
    bci["Slope",c("LCI","UCI")] <- CI.b1
    return(bci)
}
