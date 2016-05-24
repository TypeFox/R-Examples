# 2014-09-01 CJS conversion to jags
# 2012-08-30 CJS fixed problem with NAs in any() and all() in error checking
# 2011-06-13 CJS added p-values to results
# 2010-11-25 CJS pretty printing of final population estimates
# 2010-09-06 CJS forced input vectors to be vectors
# 2010-08-06 CJS added trace plots of logitP and logU
# 2010-08-03 CJS added version/date to final results
# 2010-03-12 CJS added n.sims etc to calling arguments so users can set
# 2009-12-01 CJS added openbugs/winbugs directory to argument list; added some basic error checking of arguments

TimeStratPetersenDiagErrorWHSteel_fit <-
  function( title="TSPDE-WHSteel", prefix="TSPDE-WHSteel-", 
           time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, sampfrac, 
           hatch.after=NULL, 
           bad.m2=c(), bad.u2.W.YoY=c(), bad.u2.W.1=c(), bad.u2.H.1=c(),
           logitP.cov=rep(1,length(n1)),
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
           tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05, 
           mu_xiP=logit(sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE)), 
           tau_xiP=1/var( logit((m2+.5)/(n1+1)), na.rm=TRUE),
           tauP.alpha=.001, tauP.beta=.001,
           run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted 
           debug=FALSE, debug2=FALSE,
           engine=c('jags',"openbugs")[1], 
           InitialSeed=ceiling(runif(1,min=0, max=if(engine=="jags"){1000000}else{14}))) {
# Fit a Time Stratified Petersen model with diagonal entries and with smoothing on U allowing for random error,
# covariates for the the capture probabilities, and separating the wild vs hatchery fish for STEELHEAD releases
# The steelhead are nice because 100% of hatchery fish are adipose fin clipped and no wild fish are adipose fin clipped
# The "diagonal entries" implies that no marked fish are recaptured outside the (time) stratum of release
#
   version <- '2014-09-01'
   options(width=200)

# Input parameters are
#    title - title for the analysis
#    prefix - prefix used for files created with the analysis results
#             this should be in standard Window's format, eg. JC-2002-ST-TSPDE
#             to which is appended various suffixes for plots etc
#    time   - vector of stratum numbers. For example, 9:38 would indicate that the
#             Trinity River system sampled weeks 9 to 38. If some values are omitted
#             e.g. time=10 not present, this indicates sampling did not take place this
#             week. The data are expanded and interpolation for the missing week takes place
#    n1, m2   - the input data consisting of fish marked and released and then recaptured.
#                 The n1 and m2 are used to calibrate the trap
#    u2.W.YoY - number of wild YoY fish (no clips)
#    u2.W.1   - number of wild age 1+ fish (no clips)
#    u2.H.1   - number of hatchery age 1+ fish (ad fin clipped). 100% of hatchery production is fin clipped
#    sampfrac - sampling fraction to adjust for how many days of the week was the trap operating
#              This is expressed as fraction i.e. 3 days out of 7 is expressed as 3/7=.42 etc.
#              If the trap was operating ALL days, then the SampFrac = 1. It is possible for the sampling
#              fraction to be > 1 (e.g. a mark is used for 8 days instead of 7. The data are adjusted
#              back to a 7 day week as well.
#    hatch.after - julian week AFTER which hatchery fish are released 
#    bad.m2  - list of julian numbers where the value of m2 is suspect.
#              For example, the capture rate could be extremely low.
#              These are set to NA prior to the call to OpenBugs
#    bad.u2.W.YoY
#    bad.u2.W.1
#    bad.u2.H.1 - list of julian weeks where the values of u2.* are suspect. 
#                 These are set to NA prior to the call to OpenBugs
#    logitP.cov - matrix of covariates for logit(P). If the strata times are "missing" some values, an intercept is assumed
#               for the first element of the covariance matrix and 0 for the rest of the covariates.
#               CAUTION - this MAY not be what you want to do. It is likely best to enter ALL strata
#               if you have any covariates. The default, if not specified, is a constant (the mean logit)
#    tauU.alpha, tauU.beta   - parameters for the prior on variance in spline coefficients
#    taueU.alpha, taueU.beta - parameters for the prior on variance in log(U) around fitted spline 
#    mu_xiP, tau_xiP         - parameters for the prior on mean logit(P)'s [The intercept term]
#                              The other covariates are assigned priors of a mean of 0 and a variance of 1000
#    tauP.alpha, tauP.beta   - parameters for the prior on 1/var of residual error in logit(P)'s
#    run.prob  - percentiles of run timing wanted 
#    debug  - if TRUE, then this is a test run with very small MCMC chains run to test out the data
#             and OpenBUGS will run and stop waiting for your to exit and complete

# force input vectors to be vectors
time      <- as.vector(time)
n1        <- as.vector(n1)
m2        <- as.vector(m2)
u2.W.YoY  <- as.vector(u2.W.YoY)
u2.W.1    <- as.vector(u2.W.1)
u2.H.1    <- as.vector(u2.H.1)
sampfrac  <- as.vector(sampfrac)

#  Do some basic error checking
#  1. Check that length of n1, m2, u2, sampfrac, time all match
if(var(c(length(n1),length(m2),length(u2.W.YoY),length(u2.W.1),length(u2.H.1), length(sampfrac),length(time)))>0){
   cat("***** ERROR ***** Lengths of n1, m2, u2.W.YoY, u2.W.1, u2.H.1, sampfrac, time must all be equal. They are:",
        length(n1),length(m2),length(u2.W.YoY),length(u2.W.1),length(u2.H.1), length(sampfrac),length(time),"\n")
   return()}
if(length(logitP.cov) %% length(n1) != 0){
   cat("***** ERROR ***** Dimension of covariate vector doesn't match length of n1 etc They are:",
        length(n1),length(logitP.cov),dim(logitP.cov),"\n")
   return()}
#  2. Check that m2<= n1
if(any(m2>n1, na.rm=TRUE)){
   cat("***** ERROR ***** m2 must be <= n1. The arguments are \n n1:",n1,"\n m2:",m2,"\n")
   return()}
#  3. Elements of bad.m2, bad.u2.W.YoY, bad.u2.W.1, bad.u2.H.1, and hatch.after must belong to time
if(!all(bad.m2 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.m2 must be elements of strata identifiers. You entered \n bad.m2:",bad.m2,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.W.YoY %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.W.YoY must be elements of strata identifiers. You entered \n bad.u2.W.YoY:",bad.u2.W.YoY,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.W.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.W.1 must be elements of strata identifiers. You entered \n bad.u2.W.1:",bad.u2.W.1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(bad.u2.H.1 %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** bad.u2.H.1 must be elements of strata identifiers. You entered \n bad.u2.H.1:",bad.u2.H.1,"\n Strata identifiers are \n time:",time, "\n")
   return()}
if(!all(hatch.after %in% time, na.rm=TRUE)){
   cat("***** ERROR ***** hatch.after must be elements of strata identifiers. You entered \n hatch.after:",hatch.after,"\n Strata identifiers are \n time:",time, "\n")
   return()}

   results.filename <- paste(prefix,"-results.txt",sep="")   


   sink(results.filename)
cat(paste("Time Stratified Petersen with Diagonal recaptures, error in smoothed U, separating wild and hatchery fish, STEELHEAD ONLY - ", date()))
cat("\nVersion: ", version)

cat("\n\n", title, "Results \n\n")


cat("*** Raw data *** \n")
temp<- cbind(time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, round(sampfrac,digits=2), logitP.cov)
colnames(temp)<- c('time', 'n1','m2','u2.W.YoY', 'u2.W.1', 'u2.H.1','SampFrac', paste("logitPcov[", 1:ncol(as.matrix(logitP.cov)),"]",sep="") )
print(temp) 
cat("\n\n")
cat("Hatchery fish are released AFTER strata: ", hatch.after,"\n\n")
cat("The following strata had m2       set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.W.YoY set to missing: ", 
     if(length(bad.u2.W.YoY)>0){bad.u2.W.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.W.1   set to missing: ", 
     if(length(bad.u2.W.1)>0){bad.u2.W.1} else {" NONE"}, "\n")
cat("The following strata had u2.H.1   set to missing: ", 
     if(length(bad.u2.H.1)>0){bad.u2.H.1} else {" NONE"}, "\n")


# Pooled Petersen estimator over ALL of the data including when no releases take place, bad m2, bad.u2.* values.
cat("\n\n*** Pooled Petersen Estimates based on pooling over ALL strata adjusting for sampling fractions ***\n\n")
cat("W.YoY  Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",sum(u2.W.YoY/sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.W.YoY/sampfrac, na.rm=TRUE))
cat("W.YoY  Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=", sum(u2.W.1/sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.W.1/sampfrac, na.rm=TRUE))
cat("W.1    Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("H.1    Total n1=", sum(n1, na.rm=TRUE),";  m2=",sum(m2, na.rm=TRUE),";  u2=",sum(u2.H.1/sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(n1, na.rm=TRUE), sum(m2, na.rm=TRUE), sum(u2.H.1/sampfrac, na.rm=TRUE))
cat("H.1    Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")


# Obtain the Pooled Petersen estimator for each group excluding bad.m2, bad.u2.* values
select <- (n1>0) & (!is.na(n1)) & (!is.na(m2)) & (!is.na(u2.W.YoY)) & (!is.na(u2.W.1)) & (!is.na(u2.H.1))
cat("\n\n*** Pooled Petersen Estimate excluding bad m2, u2.* values adjusting for sampling fractions ***\n\n")
cat("The following strata are excluded because n1=0 or NA values in m2, u2.W.YoY, u2.W.1, u2.H.1 :", time[!select],"\n\n")

temp.n1       <- n1      [select]
temp.m2       <- m2      [select]
temp.u2.W.YoY <- u2.W.YoY[select]
temp.u2.W.1   <- u2.W.1  [select]
temp.u2.H.1   <- u2.H.1  [select]
temp.sampfrac <- sampfrac[select]

cat("W.YoY  Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=",sum(temp.u2.W.YoY/temp.sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.W.YoY/temp.sampfrac, na.rm=TRUE))
cat("W.YoY  Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=", sum(temp.u2.W.1/temp.sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.W.1/temp.sampfrac, na.rm=TRUE))
cat("W.1    Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("H.1    Total n1=", sum(temp.n1, na.rm=TRUE),";  m2=",sum(temp.m2, na.rm=TRUE),";  u2=",sum(temp.u2.H.1/temp.sampfrac, na.rm=TRUE),"\n")
pp <- SimplePetersen(sum(temp.n1, na.rm=TRUE), sum(temp.m2, na.rm=TRUE), sum(temp.u2.H.1/sampfrac, na.rm=TRUE))
cat("H.1    Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")


# Obtain the Pooled Petersen estimator after fixup of bad.m2, bad.u2.* values after adjusting for sampling fractions
temp.m2 <- m2
index.bad.m2 <- as.vector((1:length(time)) %*% outer(time,bad.m2,"=="))
temp.m2[index.bad.m2] <- NA

temp.u2.W.YoY      <- u2.W.YoY
index.bad.u2.W.YoY <- as.vector((1:length(time)) %*% outer(time,bad.u2.W.YoY,"=="))
temp.u2.W.YoY[index.bad.u2.W.YoY] <- NA

temp.u2.W.1      <- u2.W.1
index.bad.u2.W.1  <- as.vector((1:length(time)) %*% outer(time,bad.u2.W.1,"=="))
temp.u2.W.1[index.bad.u2.W.1] <- NA

temp.u2.HJ.1      <- u2.H.1
index.bad.u2.H.1  <- as.vector((1:length(time)) %*% outer(time,bad.u2.H.1,"=="))
temp.u2.H.1[index.bad.u2.H.1] <- NA


select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) & (!is.na(temp.u2.W.YoY) & (!is.na(temp.u2.W.1)) & (!is.na(temp.u2.H.1)) )
cat("\n\n*** Pooled Petersen Estimate after removing bad m2, u2.* values adjusting for sampling fraction  ***\n\n")
cat("The following strata had m2       set to missing: ", 
     if(length(bad.m2)>0){bad.m2} else {" NONE"}, "\n")
cat("The following strata had u2.W.YoY set to missing: ", 
     if(length(bad.u2.W.YoY)>0){bad.u2.W.YoY} else {" NONE"}, "\n")
cat("The following strata had u2.W.1   set to missing: ", 
     if(length(bad.u2.W.1)>0){bad.u2.W.1} else {" NONE"}, "\n")
cat("The following strata had u2.H.1   set to missing: ", 
     if(length(bad.u2.H.1)>0){bad.u2.H.1} else {" NONE"}, "\n")
cat("\nThe following strata are excluded because n1=0 or NA values in m2, u2.*:", time[!select],"\n\n")

temp.n1        <- n1      [select]
temp.m2        <- m2      [select]
temp.u2.W.YoY  <- u2.W.YoY[select]
temp.u2.W.1    <- u2.W.1  [select]
temp.u2.H.1    <- u2.H.1  [select]
temp.sampfrac  <- sampfrac[select]

cat("W.YoY  Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.W.YoY/temp.sampfrac),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.W.YoY/temp.sampfrac))
cat("W.YoY  Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("W.1    Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.W.1/temp.sampfrac),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.W.1/temp.sampfrac))
cat("W.1    Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")

cat("H.1   Total n1=", sum(temp.n1),";  m2=",sum(temp.m2),";  u2=",sum(temp.u2.H.1/temp.sampfrac),"\n")
pp <- SimplePetersen(sum(temp.n1), sum(temp.m2), sum(temp.u2.H.1/temp.sampfrac))
cat("H.1   Est U(total) ", format(round(pp$est),big.mark=","),"  (SE ", format(round(pp$se), big.mark=","), ")\n\n\n")



# Obtain Petersen estimator for each stratum prior to removing bad m2 values
cat("*** Stratified Petersen Estimator for each stratum PRIOR to removing bad m2 values adjusting for sampling fractions ***\n\n")
cat("W.YoY  raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.W.YoY
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.YoY  Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


cat("W.1   raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.W.1
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.1    Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("H.1   raw data\n")
temp.n1 <- n1
temp.m2 <- m2
temp.u2 <- u2.H.1
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2)
temp <- cbind(time, temp.n1, temp.m2, temp.u2, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("H.1    Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")




# Obtain Petersen estimator for each stratum after removing bad m2, u2.* values
cat("*** Stratified Petersen Estimator for each stratum AFTER removing bad m2, u2.* values ***\n\n")
temp.n1 <- n1
temp.m2 <- m2
temp.m2[index.bad.m2] <- NA

temp.u2.W.YoY <- u2.W.YoY
temp.u2.W.YoY[index.bad.u2.W.YoY] <- NA

temp.u2.W.1 <- u2.W.1
temp.u2.W.1[index.bad.u2.W.1] <- NA

temp.u2.H.1 <- u2.H.1
temp.u2.H.1[index.bad.u2.H.1] <- NA

cat("\nW.YoY raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.W.YoY/sampfrac)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.W.YoY/sampfrac, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.YoY(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.YoY Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("\nW.1   raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.W.1/sampfrac)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.W.1/sampfrac, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.W.1(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("W.1   Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")

cat("\nH.1   raw data\n")
sp <- SimplePetersen(temp.n1, temp.m2, temp.u2.H.1/sampfrac)
temp <- cbind(time, temp.n1, temp.m2, temp.u2.H.1/sampfrac, round(sp$est), round(sp$se))
colnames(temp) <- c('time', 'n1','m2','u2.H.1(adj)', 'U[i]', 'SE(U[i])')
print(temp)
cat("\n")
cat("H.1   Est U(total) ", format(round(sum(sp$est, na.rm=TRUE)),big.mark=","),
    "  (SE ", format(round(sqrt(sum(sp$se^2, na.rm=TRUE))), big.mark=","), ")\n\n\n")


# Test if pooling can be done
cat("*** Test if pooled Petersen is allowable. [Check if marked fractions are equal] ***\n\n")
select <- (n1>0) & (!is.na(n1)) & (!is.na(temp.m2)) 
temp.n1 <- n1[select]
temp.m2 <- m2[select]
test <- TestIfPool( temp.n1, temp.m2)
cat("(Large Sample) Chi-square test statistic ", test$chi$statistic," has p-value", test$chi$p.value,"\n\n")
temp <- cbind(time[select],test$chi$observed, round(test$chi$expected,1), round(test$chi$residuals^2,1))
colnames(temp) <- c('time','n1-m2','m2','E[n1-m2]','E[m2]','X2[n1-m2]','X2[m2]')
print(temp)
cat("\n Be cautious of using this test in cases of small expected values. \n\n")



# Fix up any data problems and prepare for the call.
# Notice that for strata entries that are missing any covariate values, only an intercept is added

# Expand the entries in case of missing time entries
new.n1         <- rep(0, max(time)-min(time)+1)
new.m2         <- rep(0, max(time)-min(time)+1)
new.u2.W.YoY   <- rep(0, max(time)-min(time)+1)
new.u2.W.1     <- rep(0, max(time)-min(time)+1)
new.u2.H.1     <- rep(0, max(time)-min(time)+1)
new.sampfrac   <- rep(0, max(time)-min(time)+1)
new.logitP.cov <- matrix(NA, nrow=max(time)-min(time)+1, ncol=ncol(as.matrix(logitP.cov)))
new.time       <- min(time):max(time)


new.n1  [time-min(time)+1]         <- n1
new.m2  [time-min(time)+1]         <- m2
new.m2  [bad.m2-min(time)+1]       <- NA    # wipe out strata where m2 is known to be bad

new.u2.W.YoY[time-min(time)+1]             <- u2.W.YoY
new.u2.W.YoY[bad.u2.W.YoY-min(time)+1]     <- NA    # wipe out strata where u2.W.YoY is known to be bad

new.u2.W.1  [time-min(time)+1]             <- u2.W.1
new.u2.W.1  [bad.u2.W.1-min(time)+1]       <- NA    # wipe out strata where u2.W.1 is known to be bad

new.u2.H.1  [time-min(time)+1]             <- u2.H.1
new.u2.H.1  [bad.u2.H.1-min(time)+1]       <- NA    # wipe out strata where u2.W.1 is known to be bad


new.sampfrac[time-min(time)+1]   <- sampfrac
new.logitP.cov[time-min(time)+1,]<- as.matrix(logitP.cov)
new.logitP.cov[ is.na(new.logitP.cov[,1]), 1] <- 1  # insert a 1 into first columns where not specified
new.logitP.cov[ is.na(new.logitP.cov)] <- 0         # other covariates are forced to zero not in column 1


# Check for and fix problems with the data
# If n1=m2=0, then set n1 to 1, and set m2<-NA
new.m2[new.n1==0] <- NA
new.n1[new.n1==0] <- 1

# Adjust data when a stratum has less than 100% sampling fraction to "estimate" the number
# of unmarked fish that were captured. It is not necessary to adjust the n1 and m2 values 
# as these are used ONLY to estimate the capture efficiency. 
# In reality, there should be a slight adjustment
# to the precision to account for this change, but this is not done.
# Similarly, if the sampling fraction is more than 1, the adjustment forces the total unmarked catch back to a single week.
new.u2.W.YoY <- round(new.u2.W.YoY/new.sampfrac)
new.u2.W.1   <- round(new.u2.W.1  /new.sampfrac)
new.u2.H.1   <- round(new.u2.H.1  /new.sampfrac)


# Adjust for strata where sampling fraction=0. On these strata
# u2.* is set to NA so that there is NO information on U2 for this stratum
new.u2.W.YoY[new.sampfrac<.001] <- NA
new.u2.W.1  [new.sampfrac<.001] <- NA
new.u2.H.1  [new.sampfrac<.001] <- NA


# Print out the revised data
hatch.indicator <- rep('   ', max(time)-min(time)+1)
hatch.indicator[hatch.after-min(time)+1]<- '***'

cat("\n\n*** Revised data *** \n")
temp<- data.frame(time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
       sampfrac=round(new.sampfrac,digits=2), new.logitP.cov=new.logitP.cov, 
       hatch.indicator=hatch.indicator)
print(temp) 
cat("\n\n")

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
   cat("\nprior to formal call to TimeStratPetersenDiagErrorWH\n")
   browser()
}


if (debug) 
   {results <- TimeStratPetersenDiagErrorWHSteel(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
            hatch.after=hatch.after-min(time)+1, 
            logitP.cov=new.logitP.cov,
            n.chains=3, n.iter=10000, n.burnin=5000, n.sims=500,   # set to low value for debugging only
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2,  engine=engine, InitialSeed=InitialSeed)
   } else #notice R syntax requires { before the else
   {results <- TimeStratPetersenDiagErrorWHSteel(title=title, prefix=prefix, 
            time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, 
            hatch.after=hatch.after-min(time)+1, 
            logitP.cov=new.logitP.cov,
            n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.sims=n.sims, 
            tauU.alpha=tauU.alpha, tauU.beta=tauU.beta, taueU.alpha=taueU.alpha, taueU.beta=taueU.beta,
            debug=debug, debug2=debug2, engine=engine, InitialSeed=InitialSeed)
  }

# Now to create the various summary tables of the results

#browser()
# A plot of the observered log(U) on the log scale, and the final mean log(U)
plot_logU <- function(title, time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, results){
#  Plot the observed and fitted logU values along with posterior limits
#  n1, m2, u2.* are the raw data
#  results is the summary table from WinBugs

   Nstrata <- length(n1)

   u2.W.YoY[is.na(u2.W.YoY)] <- 1  # in case of missing values
   u2.W.1  [is.na(u2.W.1)]   <- 1  # in case of missing values
   u2.H.1  [is.na(u2.H.1)]   <- 1  # in case of missing values

   avg.P <- sum(m2,na.rm=TRUE)/sum(n1, na.rM=TRUE)
   Uguess.W.YoY <- pmax((u2.W.YoY+1)*(n1+2)/(m2+1), u2.W.YoY/avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
   Uguess.W.1   <- pmax((u2.W.1  +1)*(n1+2)/(m2+1), u2.W.1  /avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
   Uguess.H.1   <- pmax((u2.H.1  +1)*(n1+2)/(m2+1), u2.H.1  /avg.P, na.rm=TRUE)  # try and keep Uguess larger than observed values
   Uguess.H.1[1:(hatch.after-min(time))] <- 0   # no hatchery fish prior to release from hatchery

   min_logU <- min( log(c(Uguess.W.YoY,Uguess.W.1,Uguess.H.1)+1), na.rm=TRUE)
   max_logU <- max( log(c(Uguess.W.YoY,Uguess.W.1,Uguess.H.1)+1), na.rm=TRUE)

   # which rows contain the etaU.W[xx] and etaU.H[xx] ?
   results.row.names <- rownames(results$summary)
   etaU.row.index.W.YoY    <- grep("etaU.W.YoY", results.row.names)
   etaU.W.YoY<- results$summary[etaU.row.index.W.YoY,]
   etaU.row.index.W.1      <- grep("etaU.W.1",   results.row.names)
   etaU.W.1  <- results$summary[etaU.row.index.W.1,]
   etaU.row.index.H.1      <- grep("etaU.H.1",   results.row.names)
   etaU.H.1  <- results$summary[etaU.row.index.H.1,]
 
   min_logU <- min( c(min_logU, etaU.W.YoY[,"mean"], etaU.W.1[,"mean"], etaU.H.1[,"mean"]), na.rm=TRUE)
   max_logU <- max( c(min_logU, etaU.W.YoY[,"mean"], etaU.W.1[,"mean"], etaU.H.1[,"mean"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.YoY[,"2.5%"], etaU.W.1[,"2.5%"], etaU.H.1[,"2.5%"]), na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.YoY[,"2.5%"], etaU.W.1[,"2.5%"], etaU.H.1[,"2.5%"]), na.rm=TRUE)
   min_logU <- min( c(min_logU, etaU.W.YoY[,"97.5%"],etaU.W.1[,"97.5%"],etaU.H.1[,"97.5%"]),na.rm=TRUE)
   max_logU <- max( c(max_logU, etaU.W.YoY[,"97.5%"],etaU.W.1[,"97.5%"],etaU.H.1[,"97.5%"]),na.rm=TRUE)
   min_logU <- max(min_logU, -1)  # no point in plotting values less than 1 for abundance

   #browser()
   # the wild fish are estimated for ALL strata
   # the hatchery fish are estimated ONLY for those weeks after hatch.after
   time.W <- time
   time.H <- time>hatch.after
   # plot the raw log(U) values
   plot(time.W, log(Uguess.W.YoY), type="p", pch="w", cex=.5,
        main=paste(title,"\nFitted spline curve to raw U.W[i] U.H[i] with 95% credible intervals"),
        sub='Fitted/Smoothed/Raw values plotted for W.YoY(solid), W.1 (dash), and H.1(dash-2)',
        ylab='log(U[i])',
        xlab='Time Index', ylim=c(min_logU,max_logU))  # initial points on log scale.
   points(time.W,       log(Uguess.W.1),         pch="W", cex=.5)
   points(time[time.H], log(Uguess.H.1[time.H]), pch="h", cex=.5)

   # plot the mean of the etaU
   points(time.W   ,       etaU.W.YoY[,"mean"],       type="p", pch=19)  # fitted values for wild YoY fish
   points(time.W+.1,       etaU.W.1  [,"mean"],       type="p", pch=20)  # fitted values for wild 1+  fish
   points(time[time.H]+.2, etaU.H.1  [time.H,"mean"], type="p", pch=22)  # fitted values for hatchery fish
   lines(time.W,          etaU.W.YoY[,"mean"])  # add smoothed spline through points
   lines(time.W+.1,       etaU.W.1  [,"mean"],       lty=2)  # add smoothed spline through points
   lines(time[time.H]+.2, etaU.H.1  [time.H,"mean"], lty=3, lwd=2)

   # plot the 2.5 -> 97.5 posterior values
   segments(time.W,          etaU.W.YoY[      ,"2.5%"], time.W,          etaU.W.YoY[      ,"97.5%"])
   segments(time.W+.1,       etaU.W.1  [      ,"2.5%"], time.W,          etaU.W.1  [      ,"97.5%"], lty=2)
   segments(time[time.H]+.2, etaU.H.1  [time.H,"2.5%"], time[time.H]+.2, etaU.H.1  [time.H,"97.5%"], lty=3, lwd=2)

   # plot the spline curve before the error term is added.
   # extract the bU coefficients
   logUne.row.index.W.YoY <- grep("logUne.W.YoY", results.row.names)
   logUne.row.index.W.1   <- grep("logUne.W.1",   results.row.names)
   logUne.row.index.H.1   <- grep("logUne.H.1",   results.row.names)
 
   logUne.W.YoY<- results$summary[logUne.row.index.W.YoY,"mean"]   
   logUne.W.1  <- results$summary[logUne.row.index.W.1,  "mean"]
   logUne.H.1  <- results$summary[logUne.row.index.H.1,  "mean"]

   lines(time.W,       logUne.W.YoY,         lty=1)  # plot the curve
   lines(time.W,       logUne.W.1,           lty=2)  # plot the curve
   lines(time[time.H], logUne.H.1[time.H],   lty=3, lwd=2)  # plot the curve # too busy to plot the spline points as well

}


pdf(file=paste(prefix,"-logU.pdf",sep=""))
plot_logU(title=title, time=new.time, n1=new.n1, m2=new.m2, u2.W.YoY=new.u2.W.YoY, u2.W.1=new.u2.W.1, u2.H.1=new.u2.H.1, results=results)
dev.off()


logitP.plot <- plot_logitP(title=title, time=new.time, n1=new.n1, m2=new.m2, 
               u2=new.u2.W.YoY+new.u2.W.1+new.u2.H.1, logitP.cov=new.logitP.cov, results=results)
ggsave(plot=logitP.plot, filename=paste(prefix,"-logitP.pdf",sep=""), height=6, width=10, units="in")
results$plots$logitP.plot <- logitP.plot


# Look at autocorrelation function for Utot.W.YoY and Utot.W.1, U.tot.H.1
pdf(file=paste(prefix,"-Utot.W.YoY-acf.pdf",sep=""))
acf(results$sims.matrix[,"Utot.W.YoY"], main=paste(title,"\nAutocorrelation function for U total.W.YoY"))
dev.off()
pdf(file=paste(prefix,"-Utot.W.1-acf.pdf",sep=""))
acf(results$sims.matrix[,"Utot.W.1"], main=paste(title,"\nAutocorrelation function for U total.W.1"))
dev.off()
pdf(file=paste(prefix,"-Utot.H.1-acf.pdf",sep=""))
acf(results$sims.matrix[,"Utot.H.1"], main=paste(title,"\nAutocorrelation function for U total.H.1"))
dev.off()

# Look at the shape of the posterior distribution
pdf(file=paste(prefix,"-Utot.W.YoY-posterior.pdf",sep=""))
plot( x=density(as.vector(results$sims.array[,,"Utot.W.YoY"])), 
    main=paste(title,'\nPosterior density plot of U-total.W.YoY'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.W.YoY",c("2.5%","97.5%")])  # add vertical reference lines
dev.off()
pdf(file=paste(prefix,"-Utot.W.1-posterior.pdf",sep=""))
plot( x=density(as.vector(results$sims.array[,,"Utot.W.1"])), 
    main=paste(title,'\nPosterior density plot of U-total.W.1'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.W.1",c("2.5%","97.5%")])  # add vertical reference lines
dev.off()
pdf(file=paste(prefix,"-Utot.H.1-posterior.pdf",sep=""))
plot( x=density(as.vector(results$sims.array[,,"Utot.H.1"])), 
    main=paste(title,'\nPosterior density plot of U-total.H.1'),
    sub ="Vertical lines mark 2.5th and 97.5th percentile")
abline(v=results$summary["Utot.H.1",c("2.5%","97.5%")])  # add vertical reference lines
dev.off()


#save the Bayesian predictive distribution (Bayesian p-value plots)
pdf(file=paste(prefix,"-GOF.pdf",sep=""))
#browser()
discrep <-PredictivePosterior.TSPDE.WHSteel (time, new.n1, new.m2, new.u2.W.YoY, new.u2.W.1, new.u2.H.1, 
          expit(results$sims.list$logitP), 
          round(results$sims.list$U.W.YoY), 
          round(results$sims.list$U.W.1), 
          cbind(matrix(0,nrow=nrow(results$sims.list$U.H.1),ncol=(hatch.after-min(time)+1)),
                round(results$sims.list$U.H.1)), hatch.after) #don't forget that hatchery fish is 0 until hatch.after
gof <- PredictivePosteriorPlot.TSPDE.WHSteel(discrep)
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
# What was the initial seed
cat("\n\n*** Initial State for this run ***: ", results$Seed.initial,"\n")
cat("*** See help(modelSetRN) for details. ***\n")

# Global summary of results
cat("\n\n*** Summary of MCMC results *** \n\n")
print(results, digits.summary=3)

# Give an alternate computation of DIC based on the variance of the deviance
# Refer to http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/DIC-slides.pdf for derivation and why
# this alternate method may be superior to that automatically computed by WinBugs/OpenBugs

cat("\n\n*** Alternate DIC computation based on p_D = var(deviance)/2 \n")
results.row.names <- rownames(results$summary)
deviance.row.index<- grep("deviance", results.row.names)
deviance          <- results$summary[deviance.row.index,]
p.D <- deviance["sd"]^2/2
dic <- deviance["mean"]+p.D
cat("    D-bar: ", deviance["mean"],";  var(dev): ", deviance["sd"]^2,
    "; p.D: ", p.D, "; DIC: ", dic)

# Summary of population sizes. Add pretty printing
cat("\n\n\n\n*** Summary of Unmarked Population Size ***\n")
cat("W.YoY \n")
temp <- results$summary[ grep("Utot.W.YoY", rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nW.1   \n")
temp<- results$summary[ grep("Utot.W.1",   rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)

cat("\n\nH.1\n")
temp<- results$summary[ grep("Utot.H.1",   rownames(results$summary)),]
old.Rhat <- temp["Rhat"]
temp<- formatC(temp, big.mark=",", format="d")
temp["Rhat"] <- formatC(old.Rhat,digits=2,format="f",flag="#")
print(temp, quote=FALSE)


#browser()
time.H <- time>hatch.after
cat("\n\n\n\n*** Summary of Quantiles of Run Timing.Wild.YoY *** \n")
cat(    "    This is based on the sample weeks provided and the U.W.YoY[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.YoY, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n*** Summary of Quantiles of Run Timing.Wild.1 *** \n")
cat(    "    This is based on the sample weeks provided and the U.W.1[i] values \n") 
q <- RunTime(time=time, U=results$sims.list$U.W.1, prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n*** Summary of Quantiles of Run Timing.Hatchery.1 *** \n")
cat(    "    This is based on the sample weeks provided and the U.H.1[i] values \n") 
q <- RunTime(time=time[time.H], U=results$sims.list$U.H.1[,time.H], prob=run.prob)
temp <- rbind(apply(q,2,mean), apply(q,2,sd))
rownames(temp) <- c("Mean", "Sd")
print(round(temp,2))

cat("\n\n")
cat(paste("*** end of fit *** ", date()))

sink()


# add some of the raw data to the bugs object for simplicity in referencing it later
results$data <- list( time=time, n1=n1, m2=m2, u2.W.YoY=u2.W.YoY, u2.W.1=u2.W.1, u2.H.1=u2.H.1,
                      sampfrac=sampfrac, 
                      hatch.after=hatch.after, 
                      bad.m2=bad.m2, bad.u2.W.YoY=bad.u2.W.YoY, bad.u2.W.1=bad.u2.W.1, bad.u2.H.1=bad.u2.H.1, 
                      logitP.cov=logitP.cov,
                      version=version, date_run=date(),title=title)
results$gof <- gof

return(results)
} # end of function
