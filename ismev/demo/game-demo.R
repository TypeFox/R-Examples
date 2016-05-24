## Copyright (C) 2012 Marius Hofert and Valerie Chavez-Demoulin
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## Demo for *G*eneralized *A*dditive *M*odeling for *E*xtremes

## load packages
require(ismev)
setwd("/home/mhofert/R/2011-11-10_or_data/")

## define non-exported functions
lambda.fit <- ismev:::lambda.fit
lambda.predict <- ismev:::lambda.predict
GPD.fit <- ismev:::GPD.fit
GPD.predict <- ismev:::GPD.predict

## set seed
set.seed(271)

## frozen parameters
nyrs <- length(yrs <- 2003:2012) # years
ngrp <- length(grp <- c("A", "B")) # groups
n <- round(cbind(A = 2800 + 800 * sin(pi*seq(0.2, 1, length.out=nyrs)),
                 B = 2000 * seq(1, 2, length.out=nyrs))) # sample sizes for each year (row) and group (column)
rownames(n) <- yrs # set row names
B <- 32 # number of bootstrap replications (rather small here)
u <- 200 # threshold
edof <- 2 # equivalent degrees of freedom (since s(..., fx=FALSE), this is not fixed)
eps <- 0.005 # epsilon to determine convergence (speed up calculations)
niter <- 32 # maximal number of iterations (speed up calculations)
alpha <- 0.999 # confidence level
a <- 0.05 # significance level for confidence intervals
doPDF <- FALSE # whether graphics are saved as .pdf

## choose xi and beta depending on year and group
xi <- cbind(A   = seq(0.4, 0.6, length = nyrs),
            B   = seq(0.8, 0.2, length = nyrs))
sq <- seq(0.1, 0.5, length.out=nyrs)
beta <- cbind(A = exp(rev(sq)) / (1+xi[,"A"]),
              B = exp(seq(0.1, 0.5, length.out=nyrs)) / (1+xi[,"B"]))


### 1) A 'warm-up' to see how gam() and predict() work #########################

if(FALSE){

    set.seed(271)
    ## dummy loss data
    z <- data.frame(year  = rep(2010:2012, each=5),
                    group = sample(LETTERS[1:2], size=15, replace=TRUE),
                    loss  = rexp(15))
    ## fitting
    gamFit1 <- gam(loss ~ year + group - 1, data=z) # fit based on all provided covariates
    gamFit2 <- gam(loss ~ year - 1, data=z) # fit based on some of the provided covariates

    ## test: no covariates
    gamDummy <- gam(loss ~ 1, data=z) # fit without covariates
    lambda.fit(gamDummy) # evaluate fitted object
    (prd <- predict(gamDummy, newdata=data.frame(foo=rep(LETTERS[1:2], each=3)), se.fit=TRUE)) # => predict() should just repeat the values in this case
    lambda.predict(gamDummy) # predict from fitted object

    ## fitting results
    (res1 <- cbind(z, fit=gamFit1$fitted.values)) # same covariate combis (= (year, group)) => fits are equal
    (res2 <- cbind(z, fit=gamFit2$fitted.values)) # same covariate combis (= year) => fits are equal
    ## => gam()$fitted.values is always of length nrow(z)

    ## predict (use se.fit=TRUE to obtain standard errors [=> $fit, $se.fit])
    (newdata1 <- expand.grid(year=unique(z$year), group=levels(z$group))) # (year, group)
    (gamPred11 <- predict(gamFit1, newdata=newdata1)) # different for each (year, group) combi
    (gamPred21 <- predict(gamFit2, newdata=newdata1)) # equal for same years

    (newdata2  <- expand.grid(year=unique(z$year))) # just year
    (gamPred12  <- predict(gamFit1, newdata=newdata2)) # fails: 'group' not found
    ## => newdata needs at least the covariates which have been used in the
    ##    fitting of the model (clear since the fitting was done based on these
    ##    additional covariates)
    ## What happens if newdata contains more covariates?
    (gamPred13 <- predict(gamFit1, newdata=expand.grid(foo=1:3,
                                                       year=unique(z$year),
                                                       group=levels(z$group))))
    ## => values for each of 'foo' are equal
    (gamPred22 <- predict(gamFit2, newdata=newdata2)) # similar to gamPred11 (exactly all covariates are used in newdata)

    (newdata3 <- data.frame(year=seq(2010, 2012, by=0.5))) # intermediate values
    (gamPred13 <- predict(gamFit1, newdata=newdata3)) # fails: 'group' not found
    (gamPred23 <- predict(gamFit2, newdata=newdata3)) # works! => predict can 'interpolate'
    ## => predict() always returns vectors (fit, se.fit) of length nrow(newdata)

}


### 2) Functions ###############################################################

## example functions for generating losses
rGPD <- function(n, xi, beta) ((1-runif(n))^(-xi)-1)*beta/xi # sampling a GPD
loss <- function(n, xi, beta, u) u + rGPD(n, xi=xi, beta=beta) # sampling losses

##' @title Compute Value-at-Risk or Expected Shortfall
##' @param x matrix with three columns containing lambda, xi, and beta
##' @param alpha confidence level
##' @param u threshold
##' @param method either "VaR" for Value-at-Risk or "ES" for expected shortfall
##' @return Value-at-Risk or expected shortfall
risk.measure <- function(x, alpha, u, method=c("VaR", "ES")){
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot(ncol(x)==3)
    lambda <- x[,1]
    xi <- x[,2]
    beta <- x[,3]
    method <- match.arg(method)
    switch(method,
           "VaR"={
               u+(beta/xi)*(((1-alpha)/lambda)^(-xi)-1)
           },
           "ES"={
               ES <- (risk.measure(cbind(lambda, xi, beta),
                                   alpha=alpha, u=u, method="VaR")+beta-xi*u)/(1-xi)
               ## adjust to be Inf if xi > 1 (i.e., ES < 0)
               ## that's a convention, see p. 79 Coles (2001)
               if(any(xi > 1)) ES[xi > 1] <- Inf
               ES
           },
           stop("wrong method"))
}


### 3) Generate data ###########################################################

## recall how mapply operates on matrices: first through all rows within a column,
## then the next column etc.

## generate losses
lossList <- mapply(loss, n=n, xi=xi, beta=beta, u=u, SIMPLIFY=FALSE) # nyrs*ngrp list with corresponding sample size many losses

## determine corresponding years and covariates
ryr <- rep(yrs, ngrp)
yearList <- mapply(rep, ryr, n, SIMPLIFY=FALSE)
rgr <- rep(grp, each=nyrs)
grpList <- mapply(rep, rgr, n, SIMPLIFY=FALSE)

## build data frame containing time, covariates, and losses
## order: year 1, ..., year 1 (n[1,1]-often), year 10, ..., year 10 (n[nyrs,1]-often)
##        then c()'ed for each group
x <- data.frame(year  = unlist(yearList),
                group = unlist(grpList),
                loss  = unlist(lossList))


### 4) Model fitting ###########################################################

## Note: Exemplarily, we only fit some models to the data; a real-life
##       application would require to fit and compare (goodness-of-fit test)
##       several fitted models; see Chavez-Demoulin, Embrechts, Hofert (2014)
##       for such a study.


### 4.1) Loss frequency (standard gam() application) ###########################

### 4.1.1) lambda fitting and predicting #######################################

## determine a data set which contains the number of losses for *each* covariate
## combination
grid <- expand.grid(group=grp, year=yrs) # *all* variable combinations
lvls <- apply(grid, 1, paste, collapse=" ")
nm <- sapply(split(x$loss, factor(paste(x$group, x$year), levels=lvls)), length)
x.num <- data.frame(grid, num = nm, row.names = seq_len(length(lvls)))
x.num <- x.num[order(x.num$group, x.num$year),] # sort (simplifies VaR computation)
## => compare with n!

## fit lambda
modlam <- gam(num ~ group + s(year, by=group) - 1, data=x.num, family=poisson)

## compute fitted and predicted values incl. pointwise asymptotic CIs
lamFit <- lambda.fit(modlam)
lamPred <- lambda.predict(modlam, alpha=a)


### 4.1.2) Plot fitted and predicted lambda and CIs ############################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(lamPred$CI.low, x.num$num), max(lamPred$CI.up, x.num$num))

## setup
if(doPDF){
    file <- "game-demo_fit_lambda.pdf" # outfile name
    pdf(file=file, width=9, height=5) # open pdf device
}
## layout
layout.mat <- matrix(1:ngrp, ncol=ngrp, byrow=TRUE) # plot matrix layout
layout.mat <- rbind(layout.mat, c(ngrp+1, ngrp+2)) # add plot regions for x axis label
layout.mat <- cbind(c(ngrp+3, 0), layout.mat) # add plot regions for y axis label
layout(layout.mat, widths=c(0.2, 1, 1), heights=c(1, 0.2)) # layout
## layout.show(ngrp+3)
## plot settings
opar <- par(mar=rep.int(0,4), oma=rep.int(3,4)) # plot parameters
## plot
for(j in 1:ngrp){
    jgrp <- lamPred$covar$group==grp[j] # group boolean
    ## predicted lambda
    yr <- lamPred$covar$year[jgrp]
    plot(yr, lamPred$predict[jgrp], type="l",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n") # predicted lambda
    ## CIs
    lines(yr, lamPred$CI.low[jgrp], lty=2) # lower CI
    lines(yr, lamPred$CI.up[jgrp], lty=2) # upper CI
    ## fitted lambda
    yr <- lamFit$covar$year[jgrp] # possibly different from before
    points(yr, lamFit$fit[jgrp], pch=20) # fitted lambda
    ## actual generated number (note: lambda(t) ~= expected # of events in year t)
    points(yr, x.num[jgrp,"num"]) # use all years with non-missing data here
    ## group labels
    text(min(xlim)+0.05*diff(xlim), min(ylim)+0.95*diff(ylim),
         labels=grp[j], font=2)
}
## x axis label
plot.new()
text(0.5, 0.1, labels="Year")
## legend
plot.new()
legend(0.06, 0.35, lty=c(1,2,0), pch=c(20,NA,1), bty="n", horiz=TRUE,
       legend=c(expression(hat(lambda)),
                substitute(a.~"CI", list(a.=1-a)),
                "true number of losses"),
       text.width=strwidth("oooooooo"))
## y axis label
plot.new()
text(0.1, 0.5, srt=90,
     labels=substitute(hat(lambda)~~"with pointwise asymptotic two-sided "*a.*"% confidence intervals", list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values
if(doPDF) dev.off()


### 4.2) GPD parameters xi and beta ############################################

### 4.2.1) xi and beta fitting #################################################

## fit with bootstrap via gamGPDboot(); recall: beta = exp(nu) / (1 + xi)
sfile <- "game-demo.rds"
if(file.exists(sfile)){
    bootGPD <- readRDS(sfile) # read the bootstrapped object
} else {

    bootGPD <- gamGPDboot(x, B=B, threshold=u, datvar="loss",
                          xiFrhs = ~ group + s(year, k=edof+1, bs="cr", by=group) - 1,
                          nuFrhs = ~ group + s(year, k=edof+1, bs="cr", by=group) - 1,
                          niter=niter, epsxi=eps, epsnu=eps)
    saveRDS(bootGPD, file=sfile) # save the bootstrapped object
}

## compute fitted values of xi, beta and CIs (pointwise bootstrapped)
xibetaFit <- GPD.fit(bootGPD, alpha=a)

## compute predicted values
xibetaPred <- GPD.predict(bootGPD)


### 4.2.2) Plot fitted and predicted xi and CIs ################################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(xibetaFit$xi$CI.low, xi), max(xibetaFit$xi$CI.up, xi))
whiskex <- 0.3 # whisker extension

## setup
if(doPDF){
    file <- "game-demo_fit_xi.pdf" # outfile name
    pdf(file=file, width=9, height=5) # open pdf device
}
## layout
layout.mat <- matrix(1:ngrp, ncol=ngrp, byrow=TRUE) # plot matrix layout
layout.mat <- rbind(layout.mat, c(ngrp+1, ngrp+2)) # add plot regions for x axis label
layout.mat <- cbind(c(ngrp+3, 0), layout.mat) # add plot regions for y axis label
layout(layout.mat, widths=c(0.2, 1, 1), heights=c(1, 0.2)) # layout
## layout.show(ngrp+3)
## plot settings
opar <- par(mar=rep.int(0,4), oma=rep.int(3,4)) # plot parameters
## plot
for(j in seq_len(ngrp)){
    jgrp <- xibetaPred$xi$covar$group==grp[j] # group boolean in xibetaPred
    ## predicted xi
    plot(xibetaPred$xi$covar$year[jgrp], xibetaPred$xi$predict[jgrp], type="l",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n")
    ## CIs
    jgrp <- xibetaFit$xi$covar$group==grp[j] # group boolean in xibetaFit
    yr <- xibetaFit$xi$covar$year[jgrp]
    for(k in seq_along(yr)) {
        lines(c(yr[k], yr[k]),
              c(xibetaFit$xi$CI.low[jgrp][k], xibetaFit$xi$CI.up[jgrp][k]), lty=2) # vertical line
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(xibetaFit$xi$CI.low[jgrp][k], 2)) # lower whisker
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(xibetaFit$xi$CI.up[jgrp][k], 2)) # upper whisker
    }
    ## fitted xi
    points(yr, xibetaFit$xi$fit[jgrp], pch=20)
    ## actual xi
    points(yrs, xi[,j]) # plotting for all years
    ## group labels
    text(min(xlim)+0.95*diff(xlim), min(ylim)+0.95*diff(ylim),
         labels=grp[j], font=2)
}
## x axis label
plot.new()
text(0.5, 0.1, labels="Year")
## legend
plot.new()
legend(0.18, 0.35, lty=c(1,2,0), pch=c(20,NA,1), bty="n", horiz=TRUE,
       legend=c(expression(hat(xi)),
                substitute(a.~"CI", list(a.=1-a)),
                expression("true"~xi)),
       text.width=strwidth("oooooooo"))
## y axis label
plot.new()
text(0.1, 0.5, srt=90,
     labels=substitute(hat(xi)~~"with bootstrapped pointwise two-sided "*a.*"% confidence intervals", list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values
if(doPDF) dev.off()


### 4.2.3) Plot fitted and predicted beta and CIs ##############################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(xibetaFit$beta$CI.low, beta), max(xibetaFit$beta$CI.up, beta))

## setup
if(doPDF){
    file <- "game-demo_fit_beta.pdf" # outfile name
    pdf(file=file, width=9, height=5) # open pdf device
}
## layout
layout.mat <- matrix(1:ngrp, ncol=ngrp, byrow=TRUE) # plot matrix layout
layout.mat <- rbind(layout.mat, c(ngrp+1, ngrp+2)) # add plot regions for x axis label
layout.mat <- cbind(c(ngrp+3, 0), layout.mat) # add plot regions for y axis label
layout(layout.mat, widths=c(0.2, 1, 1), heights=c(1, 0.2)) # layout
## layout.show(ngrp+3)
## plot settings
opar <- par(mar=rep.int(0,4), oma=rep.int(3,4)) # plot parameters
## plot
for(j in 1:ngrp){
    jgrp <- xibetaPred$beta$covar$group==grp[j] # group boolean in xibetaPred
    ## predicted beta
    plot(xibetaPred$beta$covar$year[jgrp], xibetaPred$beta$predict[jgrp], type="l",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n")
    ## CIs
    jgrp <- xibetaFit$beta$covar$group==grp[j] # group boolean in xibetaFit
    yr <- xibetaFit$beta$covar$year[jgrp]
    for(k in seq_along(yr)) {
        lines(c(yr[k], yr[k]),
              c(xibetaFit$beta$CI.low[jgrp][k], xibetaFit$beta$CI.up[jgrp][k]), lty=2) # vertical line
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(xibetaFit$beta$CI.low[jgrp][k], 2)) # lower whisker
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(xibetaFit$beta$CI.up[jgrp][k], 2)) # upper whisker
    }
    ## fitted beta
    points(yr, xibetaFit$beta$fit[jgrp], pch=20)
    ## actual beta
    points(yrs, beta[,j]) # plotting for all years
    ## group labels
    text(min(xlim)+0.05*diff(xlim), min(ylim)+0.95*diff(ylim),
         labels=grp[j], font=2)
}
## x axis label
plot.new()
text(0.5, 0.1, labels="Year")
## legend
plot.new()
legend(0.18, 0.35, lty=c(1,2,0), pch=c(20,NA,1), bty="n", horiz=TRUE,
       legend=c(expression(hat(beta)),
                substitute(a.~"CI", list(a.=1-a)),
                expression("true"~beta)),
       text.width=strwidth("oooooooo"))
## y axis label
plot.new()
text(0.1, 0.5, srt=90,
     labels=substitute(hat(beta)~~"with bootstrapped pointwise two-sided "*a.*"% confidence intervals", list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values
if(doPDF) dev.off()


### 5) VaR #####################################################################

### 5.1) Fit and predict VaR and compute CIs ###################################

### fit ########################################################################

## determine how the object should look like
## lamFit, xibetaFit => fitted VaR (= function in estimators of lambda, xi, beta)
## depends on bl and year. There are 20 (group, year) combinations. However,
## for some combinations we don't have losses => only 18 combinations are available
## => use them for computing the fitted VaR.
avail <- x.num$num > 0 # boolean indicating (in x.num) which of the (group, year) combinations are available
covar <- x.num[avail, c("group", "year")] # covariate combinations for which losses are available; sorted lexicographically since x.num is
stopifnot(x.num[, c("group", "year")] == lamFit$covar) # => same order (due to sorting of x.num), we can thus use avail to pick out fitted lambda
## pick out lambda
lamFit. <- lamFit$fit[avail] # => 18
## pick out xi
stopifnot(covar == xibetaFit$xi$covar) # => same order [nothing to pick out]
xiFit.mat <- cbind(xibetaFit$xi$fit, xibetaFit$xi$boot) # => (18, B+1)
## pick out beta
stopifnot(covar == xibetaFit$beta$covar) # => same order [nothing to pick out]
betaFit.mat <- cbind(xibetaFit$beta$fit, xibetaFit$beta$boot) # => (18, B+1)

## compute fitted VaR for all those (B+1)-many vectors (lambda, xi, beta)
VaR <- cbind(covar, fit=sapply(1:(B+1), function(j)
                    risk.measure(cbind(lamFit., xiFit.mat[,j], betaFit.mat[,j]),
                                 alpha=alpha, u=u, method="VaR")))
VaR.boot <- subset(VaR, select=(ncol(VaR)-B):ncol(VaR)) # bootstrapped results
VaR.fit <- data.frame(covar, # covariates
                      fit    = VaR$fit.1, # fit
                      CI.low = apply(VaR.boot, 1, quantile, probs=a/2), # lower CI
                      CI.up  = apply(VaR.boot, 1, quantile, probs=1-a/2)) # upper CI

### predict ####################################################################

## lamPred, xibetaPred => predicted VaR (= function in predicted lambda, xi, beta)
## depends on bl and year. Predict on all 20 combinations; note: GPD.predict()
## does not (cannot) guarantee the same order of covariates as for lambda, for example
## => check! To guarantee the same order, one could either sort by hand (order())
##    or call GPD.predict() with a specific newdata (namely covar)
covar <- lamPred$covar
lamPred. <- lamPred$predict # => 20
## pick out xi
stopifnot(covar == xibetaPred$xi$covar) # => same order [nothing to pick out]
xiPred. <- xibetaPred$xi$predict # predicted beta's
## pick out beta
stopifnot(covar == xibetaPred$beta$covar) # => same order [nothing to pick out]
betaPred. <- xibetaPred$beta$predict # predicted beta's

## compute predicted VaR
VaR.pred <- data.frame(covar, # covariates
                       predict = risk.measure(cbind(lamPred., xiPred., betaPred.),
                                              alpha=alpha, u=u, method="VaR")) # predict

## basic sanity check
stopifnot(VaR.fit$CI.low < VaR.fit$fit & VaR.fit$fit < VaR.fit$CI.up)


### 5.2) VaR plot ##############################################################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(VaR.fit$CI.low), max(VaR.fit$CI.up))

## setup
if(doPDF){
    file <- "game-demo_fit_VaR.pdf" # outfile name
    pdf(file=file, width=9, height=5) # open pdf device
}
## layout
layout.mat <- matrix(1:ngrp, ncol=ngrp, byrow=TRUE) # plot matrix layout
layout.mat <- rbind(layout.mat, c(ngrp+1, ngrp+2)) # add plot regions for x axis label
layout.mat <- cbind(c(ngrp+3, 0), layout.mat) # add plot regions for y axis label
layout(layout.mat, widths=c(0.2, 1, 1), heights=c(1, 0.2)) # layout
## layout.show(ngrp+3)
## plot settings
opar <- par(mar=rep.int(0,4), oma=rep.int(3,4)) # plot parameters
## plot
for(j in 1:ngrp){
    jgrp <- VaR.pred$group==grp[j] # group boolean in VaR.pred
    ## predicted VaR
    plot(VaR.pred$year[jgrp], VaR.pred$predict[jgrp], type="l", log="y",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n")
    ## CIs
    jgrp <- VaR.fit$group==grp[j] # group boolean in VaR.fit
    yr <- VaR.fit$year[jgrp]
    for(k in seq_along(yr)) {
        lines(c(yr[k], yr[k]),
              c(VaR.fit$CI.low[jgrp][k], VaR.fit$CI.up[jgrp][k]), lty=2) # vertical line
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(VaR.fit$CI.low[jgrp][k], 2)) # lower whisker
        lines(c(yr[k]-whiskex, yr[k]+whiskex), rep(VaR.fit$CI.up[jgrp][k], 2)) # upper whisker
    }
    ## fitted VaR
    points(yr, VaR.fit$fit[jgrp], pch=20)
    ## actual VaR
    VaR.true <- risk.measure(cbind(lambda=n[,j], xi=xi[,j], beta=beta[,j]),
                             alpha=alpha, u=u)
    points(yrs, VaR.true)
    ## group labels
    text(min(xlim)+0.9*diff(xlim), min(ylim)+0.6*diff(ylim),
         labels=grp[j], font=2)
}
## x axis label
plot.new()
text(0.5, 0.1, labels="Year")
## legend
plot.new()
legend(0.1, 0.35, lty=c(1,2,0), pch=c(20,NA,1), bty="n", horiz=TRUE,
       legend=as.expression(c(substitute(widehat(VaR)[a.], list(a.=alpha)),
                substitute(a.~"CI", list(a.=1-a)),
                substitute("true"~VaR[a.], list(a.=alpha)))),
       text.width=strwidth("oooooooooo"))
## y axis label
plot.new()
text(0.1, 0.5, srt=90,
     labels=substitute(widehat(VaR)[a.]~~"with bootstrapped ptw. two-sided "*b*"% confidence intervals", list(a.=alpha, b=1-a)))
## finalize
par(opar) # reset plot parameters to their old values
if(doPDF) dev.off()


