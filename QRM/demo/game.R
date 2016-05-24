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
require(QRM)
require(mgcv)

## set seed
set.seed(271)

## frozen parameters
nyrs <- length(yrs <- 2004:2013) # years
ngrp <- length(grp <- c("A", "B")) # groups
n <- round(cbind(A = 2800 + 800 * sin(pi*seq(0.2, 1, length.out=nyrs)),
                 B = 2000 * seq(1, 2, length.out=nyrs))) # sample sizes for each year (row) and group (column)
rownames(n) <- yrs # set row names
B <- 32 # number of bootstrap replications (rather small here)
u <- 200 # threshold
edof <- 3 # degrees of freedom
## => note: by *not* specifying "fx=TRUE, k=edof+1, bs="cr" for fitting lambda, the fits are slightly better
eps <- 0.005 # epsilon to determine convergence (speed up calculations)
niter <- 32 # maximal number of iterations (speed up calculations)
alpha <- 0.999 # confidence level
a <- 0.05 # significance level for confidence intervals

## choose xi and beta depending on year and group
xi <- cbind(A   = seq(0.4, 0.6, length = nyrs),
            B   = seq(0.8, 0.2, length = nyrs))
sq <- seq(0.1, 0.5, length.out=nyrs)
beta <- cbind(A = exp(rev(sq)) / (1+xi[,"A"]),
              B = exp(seq(0.1, 0.5, length.out=nyrs)) / (1+xi[,"B"]))


### 1) A 'warm-up' to see how gam() and predict() work #########################

if(FALSE){

    require(mgcv)
    set.seed(271)
    ## dummy loss data
    z <- data.frame(year  = rep(2011:2013, each=5),
                    group = sample(LETTERS[1:2], size=15, replace=TRUE),
                    loss  = rexp(15))
    ## fitting
    ## note: ?gam says:
    ##       1) "The GAM penalized likelihood maximization problem is
    ##           solved by Penalized Iteratively Reweighted Least Squares (P-IRLS)
    ##           (see e.g. Wood 2000)"
    ##       2) Details of the default underlying fitting methods are given in
    ##          Wood (2011 and 2004)
    ##       3) gam() calls mgcv:::estimate.gam() which calls mgcv:::gam.fit()
    ##          => 'G' (first arg of mgcv:::gam.fit) contains all info
    ##          => G$mf is not used; rather G$X (defined as X) => groups are treated as 0-1
    ##             variables
    if(FALSE)
        debug(mgcv:::gam.fit)
    gamFit1 <- gam(loss ~ year + group - 1, data=z) # fit based on all provided covariates
    if(FALSE)
        undebug(mgcv:::gam.fit)
    gamFit2 <- gam(loss ~ year - 1, data=z) # fit based on some of the provided covariates

    ## test: no covariates
    gamDummy <- gam(loss ~ 1, data=z) # fit without covariates
    get.gam.fit(gamDummy) # evaluate fitted object
    (prd <- predict(gamDummy, newdata=data.frame(foo=rep(LETTERS[1:2], each=3)), se.fit=TRUE)) # => predict() should just repeat the values in this case
    gam.predict(gamDummy, value="lambda") # predict from fitted object

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

    (newdata3 <- data.frame(year=seq(2011, 2013, by=0.5))) # intermediate values
    (gamPred13 <- predict(gamFit1, newdata=newdata3)) # fails: 'group' not found
    (gamPred23 <- predict(gamFit2, newdata=newdata3)) # works! => predict can 'interpolate'
    ## => predict() always returns vectors (fit, se.fit) of length nrow(newdata)

}


### 2) Functions ###############################################################

## example functions for generating losses
loss <- function(n, xi, beta, u)
    u + rGPD(n, xi=xi, beta=beta) # sampling losses


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

## remove values for one covariate combination (to make it more interesting)
rm.y <- 2006
rm.g <- "A"
x <- x[!(x$year==rm.y & x$group==rm.g),]


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
x.num <- x.num[order(x.num$group, x.num$year),] # sort; simplifies VaR computation
## => compare with n!

## fit lambda
x.num. <- x.num[x.num$num>0,]
## modlam <- gam(num~group+s(year, fx=TRUE, k=edof+1, bs="cr")-1, # => bad fit
##               data=x.num., family=poisson)
modlam <- gam(num~group+s(year, fx=TRUE, k=edof+1, bs="cr", by=group)-1, # => fine (interaction)
              data=x.num., family=poisson)

## compute fitted and predicted values incl. pointwise asymptotic CIs
lamFit <- get.gam.fit(modlam)
lamPred <- gam.predict(modlam, alpha=a, value="lambda")


### 4.1.2) Plot fitted and predicted lambda and CIs ############################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(lamPred$CI.low, x.num$num), max(lamPred$CI.up, x.num$num))

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
    ## predicted lambda
    jgrp <- lamPred$covar$group==grp[j] # group boolean
    yr <- lamPred$covar$year[jgrp]
    plot(yr, lamPred$predict[jgrp], type="l",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n") # predicted lambda
    ## CIs
    lines(yr, lamPred$CI.low[jgrp], lty=2) # lower CI
    lines(yr, lamPred$CI.up[jgrp], lty=2) # upper CI
    ## fitted lambda
    jgrp <- lamFit$covar$group==grp[j] # possibly different from before (if there are covariate combinations without losses...)
    yr <- lamFit$covar$year[jgrp] # possibly different from before
    points(yr, lamFit$fit[jgrp], pch=20) # fitted lambda
    ## actual generated number (note: lambda(t) ~= expected # of events in year t)
    jgrp <- x.num$group==grp[j] # possibly different from before
    yr <- x.num[jgrp, "year"] # possibly different from before
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
     labels=substitute(hat(lambda)~~"with pointwise asymptotic two-sided "*a.*"% confidence intervals",
     list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values


### 4.1.3) rho fitting and predicting ##########################################

## determine a data set which contains the rate of exceedances for *each* covariate
## combination
rate.exc <- sapply(split(x$loss, factor(paste(x$group, x$year), levels=lvls)),
                   function(z) sum(z > u)/length(z)) # rates of exceedances
x.rate <- data.frame(grid, rate.exc = rate.exc, row.names = seq_len(length(lvls)))
x.rate <- x.rate[order(x.rate$group, x.rate$year),] # sort; simplifies VaR computation

## fit various models for rho and choose the 'best'
rho1 <- gam(rate.exc~1, data=x.rate, link=logit) # fit logistic regression (gam) model
rhoGrp <- gam(rate.exc~group-1, data=x.rate, link=logit) # fit gam model
1-pchisq(as.numeric(-2*(logLik(rho1)-logLik(rhoGrp))), df=1) # => group is significant

rhoYGrp <- gam(rate.exc~group+year-1, data=x.rate, link=logit) # fit gam model
1-pchisq(as.numeric(-2*(logLik(rhoGrp)-logLik(rhoYGrp))), df=1) # => year is not significant

## model for rho
modrho <- gam(rate.exc~group-1, data=x.rate, link=logit)

## compute fitted and predicted rates incl. pointwise asymptotic CIs
rhoFit <- get.gam.fit(modrho)
rhoPred <- gam.predict(modrho, alpha=a, value="rho")


### 4.2) GPD parameters xi and beta ############################################

### 4.2.1) xi and beta fitting #################################################

## fit with bootstrap via gamGPDboot(); recall: beta = exp(nu) / (1 + xi)
sfile <- "game.rds"
if(file.exists(sfile)){
    bootGPD <- readRDS(sfile) # read the bootstrapped object
} else {
    ## note: - see ?s -> by: a replicate of the smooth is produced for each factor level
    ##       - this takes some minutes... get a coffee
    ##       - the result object will be stored in 'game.rds' in your working directory (~ 800MB!)
    bootGPD <- gamGPDboot(x, B=B, threshold=u, datvar="loss",
                          xiFrhs = ~ group+s(year, fx=TRUE, k=edof+1, bs="cr", by=group)-1, # interaction
                          nuFrhs = ~ group+s(year, fx=TRUE, k=edof+1, bs="cr", by=group)-1, # interaction
                          niter=niter, eps.xi=eps, eps.nu=eps)
    saveRDS(bootGPD, file=sfile) # save the bootstrapped object (takes ~ 5min!!!)
}

## compute fitted values of xi, beta and CIs (pointwise bootstrapped)
xibetaFit <- get.GPD.fit(bootGPD, alpha=a) # takes several s

## compute predicted values
xibetaPred <- GPD.predict(bootGPD)


### 4.2.2) Plot fitted and predicted xi and CIs ################################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(xibetaFit$xi$CI.low, xi), max(xibetaFit$xi$CI.up, xi))
whiskex <- 0.3 # whisker extension

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
    ## predicted xi
    jgrp <- xibetaPred$xi$covar$group==grp[j] # group boolean in xibetaPred
    yr <- xibetaPred$xi$covar$year[jgrp]
    plot(yr, xibetaPred$xi$predict[jgrp], type="l",
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
    if(grp[j]==rm.g) points(yrs[yrs!=rm.y], xi[yrs!=rm.y, j]) else points(yrs, xi[,j])
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
     labels=substitute(hat(xi)~~"with bootstrapped pointwise two-sided "*a.*"% confidence intervals",
     list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values


### 4.2.3) Plot fitted and predicted beta and CIs ##############################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(xibetaFit$beta$CI.low, beta), max(xibetaFit$beta$CI.up, beta))

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
    ## predicted beta
    jgrp <- xibetaPred$beta$covar$group==grp[j] # group boolean in xibetaPred
    yr <- xibetaPred$beta$covar$year[jgrp]
    plot(yr, xibetaPred$beta$predict[jgrp], type="l",
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
    if(grp[j]==rm.g) points(yrs[yrs!=rm.y], beta[yrs!=rm.y, j]) else points(yrs, beta[,j])
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
     labels=substitute(hat(beta)~~"with bootstrapped pointwise two-sided "*a.*"% confidence intervals",
     list(a.=1-a)))
## finalize
par(opar) # reset plot parameters to their old values


### 5) VaR #####################################################################

### 5.1) Fit and predict VaR and compute CIs ###################################

### fit ########################################################################

## rhoFit, xibetaFit => fitted VaR (= function in estimators of rho, xi, beta)
## depends on bl and year. There are 19 available (group, year) combinations.

## note: rhoFit$covar only depends on 'group' => expand to xibetaFit$xi$covar
covar <- xibetaFit$xi$covar # => 19 (group, year) combinations
rhoFit$fit <- rhoFit$fit[match(covar$group, rhoFit$covar)]
rhoFit$covar <- covar
## pick out rho (all 1s here)
rhoFit. <- rhoFit$fit # length 19
## pick out xi
xiFit.mat <- cbind(xibetaFit$xi$fit, xibetaFit$xi$boot) # => (19, B+1)
## pick out beta
stopifnot(xibetaFit$beta$covar == covar) # => same order
betaFit.mat <- cbind(xibetaFit$beta$fit, xibetaFit$beta$boot) # => (19, B+1)

## compute fitted VaR for all those (B+1)-many vectors (rho, xi, beta)
VaR <- cbind(covar, fit=sapply(1:(B+1), function(j)
    risk.measure(cbind(rhoFit., xiFit.mat[,j], betaFit.mat[,j]),
                 alpha=alpha, u=u, method="VaR")
))
VaR.boot <- subset(VaR, select=(ncol(VaR)-B):ncol(VaR)) # bootstrapped results
VaR.fit <- data.frame(covar, # covariates
                      fit    = VaR[,"fit.1"], # fit (3rd col; first two contain covariates)
                      CI.low = apply(VaR.boot, 1, quantile, probs=a/2), # lower CI
                      CI.up  = apply(VaR.boot, 1, quantile, probs=1-a/2)) # upper CI

### predict ####################################################################

## rhoPred, xibetaPred => predicted VaR. Predict on all 20 combinations;
## note: GPD.predict() does not (cannot) guarantee the same order of covariates
##       as for rho => check!
##       To guarantee the same order, one could either sort by hand (order())
##       or call GPD.predict() with a specific newdata (namely covar)

## note: rhoPred$covar only depends on 'group' => expand to xibetaPred$xi$covar
covar <- xibetaPred$xi$covar # => all 20 (group, year) combinations
idx <- match(covar$group, rhoPred$covar$group)
rhoPred$predict <- rhoPred$predict[idx]
rhoPred$CI.low <- rhoPred$CI.low[idx]
rhoPred$CI.up <- rhoPred$CI.up[idx]
rhoPred$covar <- covar
## pick out rho (all 1s here)
rhoPred. <- rhoPred$predict # length 20
## pick out xi
xiPred. <- xibetaPred$xi$predict # predicted beta's
## pick out beta
stopifnot(xibetaPred$beta$covar == covar) # => same order
betaPred. <- xibetaPred$beta$predict # predicted beta's

## compute predicted VaR
VaR.pred <- data.frame(covar, # covariates
                       predict = risk.measure(cbind(rhoPred., xiPred., betaPred.),
                                              alpha=alpha, u=u, method="VaR")) # predict

## basic sanity check
stopifnot(VaR.fit$CI.low < VaR.fit$fit & VaR.fit$fit < VaR.fit$CI.up)


### 5.2) VaR plot ##############################################################

## plot preliminaries
xlim <- range(yrs)
ylim <- c(min(VaR.fit$CI.low), max(VaR.fit$CI.up))

## layout
doPDF <- TRUE
if(doPDF) pdf("demo_VaR_0.999.pdf", width=9, height=5)
layout.mat <- matrix(1:ngrp, ncol=ngrp, byrow=TRUE) # plot matrix layout
layout.mat <- rbind(layout.mat, c(ngrp+1, ngrp+2)) # add plot regions for x axis label
layout.mat <- cbind(c(ngrp+3, 0), layout.mat) # add plot regions for y axis label
layout(layout.mat, widths=c(0.2, 1, 1), heights=c(1, 0.2)) # layout
## layout.show(ngrp+3)
## plot settings
opar <- par(mar=rep.int(0,4), oma=rep.int(3,4)) # plot parameters
## plot
for(j in 1:ngrp){
    ## predicted VaR
    jgrp <- VaR.pred$group==grp[j] # group boolean in VaR.pred
    yr <- VaR.pred$year[jgrp]
    plot(yr, VaR.pred$predict[jgrp], type="l", log="y",
         xlim=xlim, ylim=ylim, yaxt=if(j%%2==1) "s" else "n")
    ## if(j%%2==1) eaxis(2, at=10^1:5) # => requires sfsmisc
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
    VaR.true <- risk.measure(cbind(rho=x.rate[x.rate$group==grp[j], "rate.exc"], # non-parametric estimate
                                   xi=xi[,j], beta=beta[,j]), alpha=alpha, u=u)
    points(yrs, VaR.true)
    ## group labels
    text(min(xlim)+0.9*diff(xlim), min(ylim)+0.9*diff(ylim),
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
     labels=substitute(widehat(VaR)[a.]~~"with bootstrapped ptw. two-sided "*b*"% confidence intervals",
     list(a.=alpha, b=1-a)))
## finalize
par(opar) # reset plot parameters to their old values
if(doPDF) dev.off()
