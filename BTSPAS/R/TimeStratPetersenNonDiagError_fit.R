# TO DO
#    - update documentation files
#       list of -> vector of
#       add logitP.fixed and logitP.fixed.values to parameters
#    - bayesian predictive posterior plots (the Bayesian p-values)

# 2014-09-01 CJS converstion to jags
# 2012-08-30 CJS fixed NAs problem in any() and all() in error checking
# 2011-06-13 CJS added p-values to results
# 2010-02-21 CJS changed u2 to new.u2 in expanded.m2 section
# 2010-11-25 CJS pretty printing of final population estimates
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS produced traceplot
# 2010-08-04 CJS added version/date to final result
# 2010-03-12 CJS added n.chains etc to calling arguments
# 2010-03-03 CJS allowed the user for fix some logitPs to account for time when no sampling (or other events
#                with known probability of capture take place
# 2009-12-08 CJS added some basic error checking on arguments
# 2009-12-07 CJS added bad.n1, bad.u2 arguments and fixups
# 2009-12-01 CJS added open/winbugs directory to argument list


TimeStratPetersenNonDiagError_fit <-
  function( title="TSPNDE",
           prefix="TSPNDE-", time, n1, m2, u2,
           sampfrac, jump.after=NULL,
           bad.n1=c(), bad.m2=c(), bad.u2=c(),
           logitP.cov=rep(1,length(u2)),
           logitP.fixed=NULL,
           logitP.fixed.values=NULL, n.chains=3,
           n.iter=200000, n.burnin=100000,
           n.sims=2000, tauU.alpha=1,
           tauU.beta=.05, taueU.alpha=1,
           taueU.beta=.05,
           mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)),
           tau_xiP=.6666, # need a better initial value for variation in catchability
           tauP.alpha=.001,
           tauP.beta=.001,
           run.prob=seq(0,1,.1), # what percentiles of run timing are wanted
           debug=FALSE, debug2=FALSE,
           engine=c('jags',"openbugs")[1],
           InitialSeed=ceiling(runif(1,min=0, max=if(engine=="jags"){1000000}else{14}))) {
# Fit a Time Stratified Petersen model with NON-diagonal entries and with smoothing on U allowing for random error
# This is the classical stratified Petersen model where the recoveries can take place for this and multiple
# strata later
#
    version <- '2014-09-01'
    options(width=200)

# Input parameters are
#    title  - title for the analysis (character string)
#    prefix - prefix used for files created with the analysis results
#             this should be in standard Window's format, eg. JC-2002-ST-TSPNDE
#             to which is appended various suffixes for plots etc (character string)
#    time   - vector of stratum numbers. For example, 9:38 would indicate that the
#             Trinity River system sampled weeks 9 to 38.
#             These values are never used in the analysis and only serve as labels for the weeks and for plotting purposes.
#             They should be contiguous equally spaced values and be the same length as u2.
#    n1, m2, u2 - the input data consisting of fish marked and released, recapture, and unmarked captured
#             Note that m2 is a MATRIX. The first column are those fish recaptured in the stratum of release
#             and the subsequent columns are those recoveries in later strata.
#             This is expanded to the full matrix [i,j] for released released in stratum i and recovered in stratum j
#             The vector u2 should be long enough to account for any fish that are recaptured later on
#             from releases late in the season. The bottom right diagonal of m2 may be all zeros - that is ok
#             Notice that length(u2) can be longer than length(n1)+nrow(m2).
#    sampfrac - sampling fraction to adjust for how many days of the week was the trap operating
#              This is expressed as fraction i.e. 3 days out of 7 is expressed as 3/7=.42 etc.
#              If the trap was operating ALL days, then the SampFrac = 1. It is possible for the sampling
#              fraction to be > 1 (e.g. a mark is used for 8 days instead of 7. The data are adjusted
#              back to a 7 day week as well.  This is vector of length(u2)
#    jump.after - in some cases, a single spline is still not flexible enough to cope with rapid
#                 changes in the run curve. For example, in the Trinity River project, a larger
#                 hatchery release occurs around stratum 14. This is a vector indicating the
#                 strata AFTER which the spline curve is allowed to jump.
#                 null or vector of arbitrary length.
#    bad.n1  - vector of stratum numbers where the value of n1 is suspect.
#    bad.m2  - vector of stratum numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to  OpenBugs
#    bad.u2  - vector of stratum numbers where the value of u2 is suspect.
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    logitP.fixed - vector of time values where the logitP will be specified in advance. Typically this occurs
#                   when no sampling takes place and logitP <- logit(0) = -10
#    logitP.fixed.values - vector of fixed values (on the logit scale). Use -10 for p[i] <- 0, and 10 for p[i] <- 1
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline
#    mu_xiP, tau_xiP         - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a variance of 1000
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and OpenBUGS will run and stop waiting for your to exit and complete

# force input vectors to be vectors as needed. Note that m2 is NOT a vector!
time     <- as.vector(time)
n1       <- as.vector(n1)
u2       <- as.vector(u2)
sampfrac <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time are consistent with each other.
#  In the non-diagonal case, they don't have to match
if(length(n1)!=nrow(m2)){
   cat("***** ERROR ***** Length of n1 and number of rows of m2 must be equal. They are:",
        length(n1),nrow(u2),"\n")
   return()}
if(var(c(length(u2),length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of u2, sampfrac, time must all be equal. They are:",
        length(u2),length(sampfrac),length(time),"\n")
   return()}
if(length(logitP.cov) %% length(u2) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of u2. They are:",
        length(u2),length(logitP.cov),dim(logitP.cov),"\n")
   return()}
#  2. Check that rowsum of m2<= n1
if(any(apply(m2,1,sum, na.rm=TRUE)>n1)){
   cat("***** ERROR ***** m2[i,+] must be <= n1[i]. The arguments are \n n1:",n1,"\n m2:",m2,"\n")
   return()}
#  3. Elements of bad.m2 and jump.after must belong to time
if(!all(bad.n1 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.n1 must be elements of strata identifiers. You entered \n bad.n1:",bad.n1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.m2 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",bad.m2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2 %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2 must be elements of strata identifiers. You entered \n bad.u2:",bad.u2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(jump.after %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** jump.after must be elements of strata identifiers. You entered \n jump.after:",jump.after,"\n Strata identifiers are \n time:",time, "\n")
   return()}

#  4. check that index of logitP.fixed belong to time
if(!all(logitP.fixed %in% time,na.rm=TRUE)){
   cat("***** ERROR ***** logitP.fixed must be elements of strata identifiers. You entered \n logitP.fixed:",logitP.fixed,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(length(logitP.fixed)!=length(logitP.fixed.values)){
   cat("***** ERROR ***** Lengths of logitP.fixed and logitP.fixed.values must all be equal. They are:",
        length(logitP.fixed),length(logitP.fixed.values),"\n")
   return()}


results.filename <- paste(prefix,"-results.txt",sep="")

sink(results.filename)
cat(paste("Time Stratified Petersen with Non-Diagonal recaptures, error in smoothed U, and log-normal distribution for travel time - ", date()))
cat("\nVersion: ", version)

cat("\n\n", title, "Results \n\n")

# length(sampfrac) =length(u2)
# m2(i,+) < n1(i)
# 0 < sampfrac < 1

cat("*** Raw data *** (padded to match length of u2) \n")
jump.indicator <- rep('   ', length(u2))
jump.indicator[time %in% jump.after]<- '***'
ex.n1 <- c(n1, rep(NA, length(u2)-length(n1)))
ex.m2 <- rbind(m2,matrix(NA, nrow=length(u2)-length(n1), ncol=ncol(m2)))
temp<- data.frame(time=time, n1=ex.n1, m2=ex.m2, u2=u2, sampfrac=round(sampfrac,digits=2), logitP.cov=logitP.cov, jump=jump.indicator)
print(temp)
cat("\n\n")
cat("Jump point are after strata: ", jump.after)
if(length(jump.after)==0) cat("none - A single spline is fit")
cat("\nFixed logitP indices are:", logitP.fixed)
if(length(logitP.fixed)==0) cat("none - NO fixed values")
cat("\nFixed logitP values  are:", logitP.fixed.values)
if(length(logitP.fixed)==0) cat("none - NO fixed values")


# Obtain the Pooled Petersen estimator prior to fixup of bad.n1, bad.m2, and bad.u2 values
cat("\n\n*** Pooled Petersen Estimate prior to fixing bad n1, m2, or u2 values  CHECK - CHECK - CHECK - CHECK ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2

cat("Total n1=", sum(temp.n1,na.rm=TRUE),";  m2=",sum(temp.m2,na.rm=TRUE),";  u2=",sum(temp.u2,na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1,na.rm=TRUE), sum(temp.m2,na.rm=TRUE), sum(temp.u2,na.rm=TRUE))
cat("Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

# Obtain the Pooled Petersen estimator after removal of entries with bad.n1, m2, or u2 values
# select <- !(time %in% bad.n1 | time %in% bad.m2 | time %in% bad.u2)
select <- (temp.n1>0) & (!is.na(n1)) & (!apply(is.na(temp.m2),1,any)) & (!is.na(temp.u2[1:length(n1)]))
cat("\n\n*** Pooled Petersen Estimate after fixing bad m2 values  CHECK - CHECK - CHECK - CHECK ***\n\n")
cat("The following strata were excluded:",
     if(length(time[!select])>0){time[!select]} else {" NONE"}, "\n")

#temp.n1 <- n1[select]
#temp.m2 <- m2[select]
#temp.u2 <- u2[select]

cat("Total n1=", sum(temp.n1,na.rm=TRUE),";  m2=",sum(temp.m2,na.rm=TRUE),";  u2=",sum(temp.u2, na.rm=TRUE),"\n\n")
pp <- SimplePetersen(sum(temp.n1,na.rm=TRUE), sum(temp.m2,na.rm=TRUE), sum(temp.u2,na.rm=TRUE))
cat("Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")


######## This needs more thoughs #########
# Obtain Petersen estimator for each stratum prior to removing bad m2 values
#cat("*** Stratified Petersen Estimator for each stratum PRIOR to removing bad m2 values ***\n\n")
#temp.n1 <- n1
#temp.m2 <- m2
#temp.u2 <- u2
#sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
#temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
#colnames(temp) <- c('time', 'n1','m2','u2', 'U[i]', 'SE(U[i])')
#print(temp)
#cat("\n")
#cat("Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
#   "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


# Obtain Petersen estimator for each stratum after removing bad m2 values
#cat("*** Stratified Petersen Estimator for each stratum AFTER removing bad m2 values ***\n\n")
#temp.n1 <- n1
#temp.m2 <- m2
#temp.m2[index.bad.m2] <- NA
#temp.u2 <- u2
#sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
#temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
#colnames(temp) <- c('time', 'n1','m2','u2', 'U[i]', 'SE(U[i])')
#print(temp)
#cat("\n")
#cat("Est U(total) ", format(round(sum(sp$est)),big.mark=","),
#    "  (SE ", format(round(sqrt(sum(sp$se^2))), big.mark=","), ")\n\n\n")




####### This needs more thought #############
# Test if pooling can be done
#cat("*** Test if pooled Petersen is allowable. [Check if marked fractions are equal] ***\n\n")
#select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2))
#temp.n1 <- n1[select]
#temp.m2 <- m2[select]
#test <- TestIfPool( temp.n1, temp.m2)
#cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
#temp <- cbind(time[select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
#colnames(temp) <- c('time','n1-m2','m2','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
#print(temp)
#cat("\n Be cautious of using this test in cases of small expected values. \n\n")


# Adjust the data for the explicity bad values or other problems
new.time <- time
new.n1   <- n1
new.m2   <- m2
new.u2   <- u2
new.sampfrac <- sampfrac
new.logitP.cov <- logitP.cov

# If n1=m2=0, then set n1 to 1, and set m2<-NA as winbugs cannot deal with n1=0 and m2=0
new.m2[new.n1==0,] <- NA
new.n1[new.n1==0 ] <- 1

########## This needs more thought ##########
# Adjust data when a stratum has less than 100% sampling fraction to "estimate" the number
# of unmarked fish that were captured. It is not necessary to adjust the n1 and m2 values
# as these are used ONLY to estimate the capture efficiency.
# In reality, there should be a slight adjustment
# to the precision to account for this change, but this is not done.
# Similarly, if the sampling fraction is more than 1, the adjustment is made back to a standard week.
#new.u2 <- round(new.u2/new.sampfrac)

# Adjust for strata where sampling fraction is very small, but the logitP has not been fixed
# to zero. On these strata
# u2 is set to NA so that there is NO information on U2 for this stratum
new.u2[new.sampfrac<.001 & new.sampfrac > 0] <- NA

# Set the bad values to missing
new.n1[time[1:length(n1)] %in% bad.n1]  <- NA
new.m2[time[1:length(n1)] %in% bad.m2,] <- NA
new.u2[time %in% bad.u2]                <- NA

# Print out the revised data
cat("\n\n*** Revised data *** \n")
jump.indicator <- rep('   ', length(u2))
jump.indicator[time %in% jump.after]<- '***'
ex.n1 <- c(new.n1, rep(NA, length(new.u2)-length(new.n1)))
ex.m2 <- rbind(new.m2,matrix(NA, nrow=length(new.u2)-length(new.n1), ncol=ncol(new.m2)))
temp<- data.frame(time=new.time, n1=ex.n1, m2=ex.m2, u2=new.u2, sampfrac=round(new.sampfrac,digits=2), logitP.cov=new.logitP.cov,
       jump.after=jump.indicator)
print(temp)
cat("\n\n")

# Need to expand the m2 matrix to the full Nstrata.rel x Nstrata.cap+1 matrix with
# the [i,j] entry for the number of fish released in stratum i and recaptured in stratum j
# The last column of the expanded m2 matrix is the number of fish released but never recaptured later.

expanded.m2 <- matrix(0, nrow=length(new.n1), ncol=length(new.u2)+1)
for(i in 1:length(new.n1)){
   expanded.m2[i,1:length(new.u2)] <- c(rep(0,i-1),new.m2[i,],rep(0,length(new.u2)))[1:length(new.u2)]
   expanded.m2[i,length(new.u2)+1] <- new.n1[i] - sum(new.m2[i,])
}

cat("*** Expanded m2 array ***\n\n")
print(expanded.m2)

# assign the logitP fixed values etc.
new.logitP.fixed <- rep(NA, length(new.u2))
new.logitP.fixed[match(logitP.fixed, time)] <- logitP.fixed.values


# Print out information on the prior distributions used
cat("\n\n*** Information on priors *** \n")
cat("   Parameters for prior on tauU (variance in spline coefficients: ", tauU.alpha, tauU.beta,
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauU.alpha/tauU.beta,2),round(sqrt(tauU.alpha/tauU.beta^2),2),"\n")
cat("   Parameters for prior on taueU (variance of log(U) about spline: ",taueU.alpha, taueU.beta,
    " which corresponds to a mean/std dev of 1/var of:",
    round(taueU.alpha/taueU.beta,2),round(sqrt(taueU.alpha/taueU.beta^2),2),"\n")
cat("   Parameters for prior on beta.logitP[1] (intercept) (mean, 1/var):", round(mu_xiP,3), round(tau_xiP,5),
    " which corresponds to a median P of ", round(expit(mu_xiP),3), "\n")
cat("   Parameters for prior on tauP (residual variance of logit(P) after adjusting for covariates: ",tauP.alpha, tauP.beta,
    " which corresponds to a mean/std dev of 1/var of:",
    round(tauP.alpha/tauP.beta,2),round(sqrt(tauP.alpha/tauP.beta^2),2),"\n")


cat("\n\nInitial seed for this run is: ",InitialSeed, "\n")

sink()

if (debug2) {
   cat("\nprior to formal call to TimeStratPetersenNonDiagError\n")
   browser()
}

if (debug)
  {results <- TimeStratPetersenNonDiagError(title=title, prefix=prefix,
                                            time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2,
                                            jump.after=(1:length(u2))[time %in% jump.after],
                                            logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
                                            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,  # set to small values for debugging only
                                            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
                                            debug=debug, debug2=debug2, engine=engine, InitialSeed=InitialSeed)
 }
else #notice R syntax requires { before the else
  {results <- TimeStratPetersenNonDiagError(title=title, prefix=prefix,
                                            time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2,
                                            jump.after=(1:length(u2))[time %in% jump.after],
                                            logitP.cov=new.logitP.cov, logitP.fixed=new.logitP.fixed,
                                            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims,
                                            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
                                            debug=debug, debug2=debug2, engine=engine, InitialSeed=InitialSeed)
 }

# Now to create the various summary tables of the results

Nstrata.rel <- length(n1)
Nstrata.cap <- ncol(expanded.m2) -1 # don't forget that last column of m2 is number of fish never seen

# A plot of the observered log(U) on the log scale, and the final mean log(U)
plot_logU <- function(title, time, n1, m2, u2, logitP.cov, results){
#  Plot the observed and fitted logU values along with posterior limits
#  n1, m2 (the expanded version), u2 are the raw data
#  results is the summary table from WinBugs

   Nstrata.rel <- length(n1)
   Nstrata.cap <- ncol(m2)-1  # remember that last column of m2 is number of fish never seen again
   Uguess <- (u2[1:Nstrata.rel]+1)*(n1+2)/(apply(m2[,1:Nstrata.cap],1,sum)+1)  # try and keep Uguess larger than observed values
   plot(time[1:Nstrata.rel], log(Uguess), type="p",
        main=paste(title,"\nFitted spline curve to raw U[i] with 95% credible intervals"),
        sub='Open/closed circles - initial and final estimates',
        ylab='log(U[i])',
        xlab='Stratum',xlim=c(min(time),max(time)))  # initial points on log scale.


   # which rows contain the etaU[xx] ?
   results.row.names <- rownames(results$summary)
   etaU.row.index    <- grep("etaU", results.row.names)
   etaU<- results$summary[etaU.row.index,]

   # plot the mean of the etaU
   points(time, etaU[,"mean"], type="p", pch=19)  # fitted values
   lines(time, etaU[,"mean"])  # add smoothed spline through points
   # plot the 2.5 -> 97.5 posterior values
   segments(time, etaU[,"2.5%"], time, etaU[,"97.5%"])

   # plot the spline curve before the error term is added.
   # extract the bU coefficients
   logUne.row.index <- grep("logUne", results.row.names)
   logUne<- results$summary[logUne.row.index,"mean"]
   points(time, logUne, type="p", pch=20)
   lines(time, logUne, lty=2)  # plot the curve
}



pdf(file=paste(prefix,"-logU.pdf",sep=""))
plot_logU(title=title, time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2, results=results)
dev.off()

logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=expanded.m2, u2=new.u2, logitP.cov=new.logitP.cov, results=results)
ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot

# Look at autocorrelation function for Ntot
pdf(file=paste(prefix,"-Utot-acf.pdf",sep=""))
acf(results$sims.matrix[,"Utot"], main=paste(title,"\nAutocorrelation function for U total"))
dev.off()

# Look at the shape of the posterior distribution
pdf(file=paste(prefix,"-Ntot-posterior.pdf",sep=""))
plot( x=density(as.vector(results$sims.array[,,"Ntot"])),
    main=paste(title,'\nPosterior density plot of N-total'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Ntot",c("2.5%","97.5%")])  # add vertical reference lines
dev.off()

pdf(file=paste(prefix,"-Utot-posterior.pdf",sep=""))
plot( x=density(as.vector(results$sims.array[,,"Utot"])),
    main=paste(title,'\nPosterior density plot of U-total'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot",c("2.5%","97.5%")])  # add vertical reference lines
dev.off()

# plot the mean log(travel times) (the muLogTT) vs release stratum number
pdf(file=paste(prefix,"-muLogTT.pdf",sep=""))
# which rows contain the muLogTT[xx] ?
results.row.names <- rownames(results$summary)
muLogTT.row.index    <- grep("muLogTT", results.row.names)
muLogTT<- results$summary[muLogTT.row.index,]
ylim <- c( min(c(muLogTT[,"mean"],muLogTT[,"2.5%"] ), na.rm=TRUE),
           max(c(muLogTT[,"mean"],muLogTT[,"97.5%"]), na.rm=TRUE))
plot (time[1:Nstrata.rel], muLogTT[,"mean"], type="p", pch=19,
      main=paste(title,"\n mean log(travel time) with 95% credible intervals"),
      xlab="Stratum", ylab="mean log(travel time)", ylim=ylim)  # fitted values
lines(time[1:Nstrata.rel],  muLogTT[,"mean"])  # join the points
segments(time[1:Nstrata.rel], muLogTT[,"2.5%"], time[1:Nstrata.rel], muLogTT[,"97.5%"])  # plot the 2.5 -> 97.5 posterior values
dev.off()

# plot the sd log(travel times) (the sdLogTT) vs release stratum number
pdf(file=paste(prefix,"-sdLogTT.pdf",sep=""))
# which rows contain the sdLogTT[xx] ?
results.row.names <- rownames(results$summary)
sdLogTT.row.index    <- grep("sdLogTT", results.row.names)
sdLogTT<- results$summary[sdLogTT.row.index,]
ylim <- c( min(c(sdLogTT[,"mean"],sdLogTT[,"2.5%"] ), na.rm=TRUE),
           max(c(sdLogTT[,"mean"],sdLogTT[,"97.5%"]), na.rm=TRUE))
plot(time[1:Nstrata.rel], sdLogTT[,"mean"], type="p", pch=19,
      main=paste(title,"\n sd log(travel time) with 95% credible intervals"),
      xlab='Stratum', ylab="sd(log travel time)", ylim=ylim)  # fitted values
lines(time[1:Nstrata.rel],  sdLogTT[,"mean"])  # join the points
segments(time[1:Nstrata.rel], sdLogTT[,"2.5%"], time[1:Nstrata.rel], sdLogTT[,"97.5%"])  # plot the 2.5 -> 97.5 posterior values
dev.off()

## save the Bayesian predictive distribution (Bayesian p-value plots)

pdf(file=paste(prefix,"-GOF.pdf",sep=""))

discrep <-PredictivePosterior.TSPNDE (new.n1, expanded.m2, new.u2,
                                      new.logitP.fixed,
                                      expit(results$sims.list$logitP),
                                      round(results$sims.list$U),
                                      results$sims.list$muLogTT,
                                      results$sims.list$sdLogTT,
				      engine)
#browser()
gof <- PredictivePosteriorPlot.TSPNDE (discrep)

dev.off()

varnames <- names(results$sims.array[1,1,])  # extract the names of the variables
# First do the trace plots of logitP
pdf(file=paste(prefix,"-trace-logitP.pdf",sep=""))
parm.names <- varnames[grep("^logitP", varnames)]
trace_plot(title=title, results=results,
    parms_to_plot=parm.names, panels=c(3,2))
dev.off()

# now for the traceplots of logU (etaU), Utot, and Ntot
pdf(file=paste(prefix,"-trace-logU.pdf",sep=""))
parm.names <- varnames[c(grep("Utot",varnames), grep("Ntot",varnames), grep("^etaU", varnames))]
trace_plot(title=title, results=results,
    parms_to_plot=parm.names, panels=c(3,2))
dev.off()

sink(results.filename, append=TRUE)

# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
print(results, digits.summary=3)

cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
results.row.names <- rownames(results$summary)
deviance.row.index<- grep("deviance", results.row.names)
deviance          <- results$summary[deviance.row.index,]
p.D <- deviance["sd"]^2/2
dic <- deviance["mean"]+p.D
cat("    D-bar: ", deviance["mean"],";  var(dev): ", deviance["sd"]^2,
    "; p.D: ", p.D, "; DIC: ", dic)

# Summary of population sizes. Add pretty printing of results.
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
temp<- results$summary[ grep("Utot", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\n*** Summary of Total Population Size *** \n")
temp<- results$summary[ grep("Ntot", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)



cat("\n\n\n\n*** Summary of Quantiles of Run Timing *** \n")
cat(    "    This is based on the sample weeks provided and the U[i] values \n")
q <- RunTime(time=time, U=results$sims.list$U, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))


cat("\n\n")
cat(paste("*** end of fit *** ", date()))

sink()


# add some of the raw data to the bugs object for simplicity in referencing it later
results$data <- list( time=time, n1=n1, m2=m2, u2=u2, sampfrac=sampfrac,
                      jump.after=jump.after,
                      bad.n1=bad.n1, bad.m2=bad.m2, bad.u2=bad.u2,
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(),title=title)
results$gof <- gof

return(results)
} # end of function

