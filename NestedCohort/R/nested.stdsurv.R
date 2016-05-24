nested.stdsurv <-
  ##
  # Compute Cox model survival and PAR, standardized for confounders, for studies nested
  # within cohorts.
  # By: Hormuzd Katki 4/10/09
  #
  #   exposures: The exposures of interest, must be factor variables
  # confounders: The confounders of interest, must be factor variables
  # Don't use '*' in exposures or confounders, use interaction()
  # IF A VARIABLE IS CATEGORICAL, IT MUST BE OF TYPE FACTOR, EVEN IF IT IS BINARY
  # DO NOT USE ANY STRATA() AND CLUSTER() AND OFFSET() STATEMENTS
  # Thus you cannot stratify the baseline hazard.
  # No left truncation, no late entry -- all must enter at 'time zero'.
  # Reason why covariate is missing cannot depend on time (no incidence-density sampling).
  # Does not account for competing risks.
  # Requires a data frame, name no variables in it as 'o.b.s.e.r.v.e.d.' or 'p.i.h.a.t.'
  # The default exposureofinterest of "" means choose the first level of the exposure.
  # Note that standardizing over nothing yields all the covariate-specific survivals.
  ##
  function(outcome,exposures,confounders, samplingmod, data,
           exposureofinterest = "", timeofinterest = Inf, cuminc = FALSE,
           plot=FALSE, plotfilename="",
           glmlink =binomial(link="logit"),
           glmcontrol = glm.control(epsilon=1e-10, maxit=10, trace=FALSE),
           coxphcontrol = coxph.control(eps=1e-10,iter.max=50),
           missvarwarn = TRUE,...)
{

# Need generalized inverses from MASS
library(MASS)

# Check if data.frame is entered
if ( is.null(data) | !is.data.frame(data) )
  stop("data is a required argument and must be a data.frame")

# Make sure exposures and confounders have no '*',or strata(),cluster(),offset() statements
illegal <- c("*","strata(","cluster(","offset(")
if ( any(unlist(lapply(illegal,grep,x=paste(exposures,confounders),fixed=TRUE))) )
  stop(paste("exposures or confounders contain some illegal commands:",
             paste(illegal,collapse=", "),
             "\n Remember: instead of using '*' use 'interaction()'"))

# Check that there are no variables in the dataset named "o.b.s.e.r.v.e.d."
if ( any("o.b.s.e.r.v.e.d." == colnames(data)) )
  stop("You are not allowed to name any variable as 'o.b.s.e.r.v.e.d.'.")

# Check that there are no variables named "p.i.h.a.t."
if ( any("p.i.h.a.t." == colnames(data)) )
  stop("You are not allowed to name any variable in the data as 'p.i.h.a.t.'.")

# Create the Cox formula by pasting together the outcome to confounders to exposures
coxformula <- paste(c(outcome,"~",confounders,"+",exposures),collapse="")

# Generate vector of which observations have a fully-observed covariate vector
# This is 1 if all exposures and confounders are observed, 0 otherwise
# Then make the formula for the sampling model
alllist <- eval(attr(terms(as.formula(coxformula)),"variables"),envir=data)
covarlist <- alllist[-1]
outcomelist <- alllist[1]
observed <- as.numeric(complete.cases(covarlist))

# Extract the Surv() object out of the model formula, keep only the observed
# First column is the lifetime, second column the event indicator
# Make sure everyone enters cohort simultaneously
if ( dim(outcomelist[[1]])[2] > 2 )
  stop("Late entry not supported: all subjects must enter cohort at 'time zero'")

# Put observed into the dataframe and create the sampling model formula
data$o.b.s.e.r.v.e.d. <- observed
if (samplingmod != "")
  samplingformula <- paste("o.b.s.e.r.v.e.d. ~",samplingmod)
else # If nothing to stratify on, then just regress on the intercept
  samplingformula <- paste("o.b.s.e.r.v.e.d. ~ 1")
  
# Make sure that exposures and confounders are factors
if ( any(sapply(covarlist,inherits,"factor")==FALSE) )
  stop("All exposures and confounders must inherit from class factor")

# If you are missing data in the samplingmod or outcome, get rid of it,
# but issue a warning to the user
temp <- cbind(outcomelist,
              as.data.frame(eval(attr(terms(as.formula(samplingformula)),"variables"),
                                 envir=data)))
missing <- apply(is.na(temp),1,any)
if ( any(missing) ) {
  if ( missvarwarn ) {
    warning(paste("You had",sum(missing),
                  "obervations with missing values in the samplingmod",
                  "or survival outcome.\n",
                  "These observations were removed.\n",
                  "There shouldn't be any missing data in the samplingmod or survival",
                  "outcome, so be careful.\n"
                  ))
  }
  data <- data[!missing,]

  # Regenerate the covarlist and who's observed:
  
  # Generate vector of which observations have a fully-observed covariate vector
  # Note that I don't need to regenerate o.b.s.e.r.v.e.d. since that's taken care of in
  # covarlist.  This is 1 if all exposures and confounders are observed, 0 otherwise
  covarlist <- eval(attr(terms(as.formula(coxformula)),"variables"),envir=data)[-1]
  observed <- as.numeric(complete.cases(covarlist))
}

# Cohort size
n <- dim(data)[1]

# I must use treatment and polynomial contrasts, so set these temporarily, then revert
# back to the user's settings.  I got this code from the Package haplo.stats()
contrasts.tmp <- c(factor = "contr.treatment", ordered = "contr.poly")
contrasts.old <- options()$contrasts
options(contrasts = contrasts.tmp)
on.exit(options(contrasts = contrasts.old))

# First fit sampling model; remember to keep the model matrix by setting x=TRUE.
# Make sure that no covariate combination has zero chance of being observed
suppressWarnings(
                 samplingmod <- glm(as.formula(samplingformula), x=TRUE, family=glmlink,
                                    na.action=na.fail, data=data, control=glmcontrol,...)
                 )

# Get the estimated pihats to use for weighting.
pihat <- data$p.i.h.a.t. <- predict(samplingmod,type='response')

# If any pihats are exactly zero, stop right here -- cannot make estimates!
if ( any(pihat==0) )
  stop(paste("According to the sampling model,",
             "some subjects have zero probability of being observed. \n",
             "Check the sampling model and try again.\n")
       )

# Check if the sampling model converged, and if it didn't, check if it's because some
# pihats are close to 1 or 0.  We expect that if some pihats are close to 1, the routine
# will never claim that it has converged, but it's fine if some pihats are close to 1.
# But if all pihats are not close to either 0 or 1 (close enough as determined by
# criticalpihat), then there's a meaningful lack of convergence that we need to warn the
# user about.
criticalpihat <- .999
if ( samplingmod$converged==FALSE )
  if ( all(pihat<criticalpihat) & all(pihat>1-criticalpihat) ) {
    warning("The sampling model didn't converge. Outputing sampling model for inspection.")
    outputsamplingmod=TRUE
  }


# Get the scores to use for the variance of the betas with 1/pihat weights.
scores <- matrix(observed-pihat,nrow=n,ncol=length(samplingmod$coeff)) * samplingmod$x

# Fit Cox model, weight by 1/pihat, don't use robust=TRUE or cluster() so that you don't
# mislead yourself into thinking that the robust variance is the correct variance.  coxmod
# <- coxph(as.formula(coxformula), weights=1/pihat,method="breslow",na.action=na.omit,x=TRUE,
# data=data,control=coxph.control(eps=1e-10,iter.max=50),...)
coxmod <- coxph(as.formula(coxformula), weights=1/pihat,method="breslow",
                na.action=na.omit, x=TRUE, y=TRUE, subset=TRUE , control=coxphcontrol, data=data)

# Get the influence function for the betahats, set to zero for everyone with missing
# covariates.
D2 <- matrix(0,nrow=n,ncol=length(coxmod$coeff))
D2[observed==1,]  <- residuals(coxmod,'dfbeta',weighted=TRUE)

# Compute the pihat subtract-off by projecting the influences onto the scores.
# Pseudoinverse is important if some combo of J*Delta is always/never observed.
gamma2 <- scores %*% ginv(t(scores)%*%scores) %*% (t(scores) %*% D2)

# This is the correct variance of betahat that accounts for the estimated weights.
D2pihat <- D2-gamma2
varbeta <- t(D2pihat)%*%(D2pihat)

# Output Cox model results, by plugging in the correct variance into coxmod's var field
# But remember the pure inverse information just in case you want it (don't call it
# naive.var, that's for robust=TRUE only, and coxph.print() will print it out also, which we
# don't want)
coxmod$wrongvar <- coxmod$var
coxmod$var <- varbeta

##
# Compute the overall Wald test
# The LRT cannot be computed since this is a pseudolikelihood
# The score test is not implemented
##

# Invalidate the loglik and score 
coxmod$loglik <- NA
coxmod$score <- NA

# Plug in the correct Wald test
# This is the usual Wald test, variance estimated at the alternative
# Remove any aliased columns, which will have NA for that coefficient
not.aliased <- which(!is.na(coxmod$coeff))
coeff <- coxmod$coeff[not.aliased]
coxmod$wald.test <- t(coeff) %*% ginv(varbeta[not.aliased,not.aliased]) %*% coeff

# Print out
summary(coxmod)



##
# End relative risk computations, now compute standardized survival and risk differences,
# and PAR
##


# Keep beta and keep observations with complete data
beta <- coxmod$coef
fulldata <- data[observed==1,]
# Extract the Surv() object out of the cox model formula
Survfulldata <- eval(attr(terms(as.formula(coxformula)),"variables"),envir=fulldata)[[1]]
lifetime <- Survfulldata[,1]
event <- Survfulldata[,2]

# Get the details of the Cox model fit and extract needed quantities Extract out the time,
# S0, and weighted-dNbar, and remember to add in time=0, S0(time=0), dNbar(time=0)=0.
# Note that the hazard of coxph.details() is the centered hazard, it's the baseline hazard
# multiplied by the exp(unweightedmeancovariate * beta) Note that coxph.detail gives an
# unintelligible warning, I don't think it matters, so I suppress it.
coxdetails  <- suppressWarnings(coxph.detail(coxmod))
time <- c(0,coxdetails$time)
hazard <- coxdetails$hazard/exp(apply(coxdetails$x,2,mean) %*% beta)

# Get weighted death process dNbar Sum up all the weights (order by event-time) for only
# the events.  Deal with tied event times by sending the reversed event times to the
# duplicated function (which picks out the times that are duplicates of times with smaller
# subscripts), then reverse back and negate.  Thus we pick out the appropriate elements of
# the cumsum.  This is Nbar, so difference to get dNbar.  Finally, add a zero at the
# beginning for time zero.
Nbar <- cumsum(coxdetails$weights[coxdetails$y[,3]==1])[!rev(duplicated(rev(coxdetails$y[coxdetails$y[,3]==1,2])))]
dNbar <- c(0,Nbar[1],diff(Nbar))
# Now back out S0 out of the hazard and dNbar.  Remember to add in the value at time zero.
S0 <- c( sum(coxdetails$weight * exp(coxdetails$x %*% beta)) , dNbar[-1] / hazard  )

# Quickly check if timeofinterest is a legal time, and if it's Inf, set to max(time)
usertimeofinterest <- timeofinterest
if (is.infinite(timeofinterest))
  timeofinterest=max(time)
else if (!is.numeric(timeofinterest) | timeofinterest<0 | timeofinterest>max(time))
  stop("Time of interest must be between 0 and end of follow-up (use Inf to mean end of follow-up)")


# This is the baseline survival estimate at last event time, which exponentiates the
# baseline cumhaz
cumhaz <- sum(hazard[time[-1]<=timeofinterest])

# Compute riskset (the index of the last failure time at risk) for each subject.
riskset <- matrix(0,nrow=dim(fulldata)[1],ncol=1)

# This below line works in R, but SPLUS doesn't have the which.max function -- which finds
# the the index of the maximum of a vector.  If there are ties, it takes the first index.
# I want the last index which is why I have to use rev to reverse the vector.
# Use the below 2 lines in SPLUS, workaround for not having the which.max function
# for (i in 1:dim(fulldata)[1]) {
#   eligible <- rev(time <= fulldata$X[i])
#   riskset[i] <- (length(time)+1)-which(eligible==max(eligible))[1]}

for (i in 1:dim(fulldata)[1]) {
  riskset[i] <- (length(time)+1)-which.max(rev(time <= lifetime[i]))
}

# Compute influence function for lambda given fixed beta (D_i^F*R_i/pi_i) at the time of
# interest.  Make sure that all those R_i=0 have D1=0 also. 
#risksetofinterest <- length(S0)
risksetofinterest <- sum(time<=timeofinterest) 
correction <- cumsum(S0^-2 * dNbar)
D1 <- matrix(0,nrow=n,ncol=1)
D1[observed==1,] <-(1/fulldata$p.i.h.a.t.)*((1/S0[riskset])*event*
                   (riskset<=risksetofinterest)
                   - correction[pmin(riskset,risksetofinterest)]*exp(coxmod$x %*% beta) )

# Compute influence function for lambda D_i^3 First compute integral of E*dLambda up to
# the time of interest. Note that if the last observed time
# is not of interest, replace the EdLambda <- line below with:
E <- coxdetails$mean
r <- risksetofinterest-1 # subtract 1 since neither E nor hazard includes time zero
EdLambda <- apply( E[1:r,] * matrix(hazard[1:r],nrow=r,ncol=dim(coxmod$x)[2]),2,sum )
D3 <- D1 - ( D2 %*% EdLambda )

# Compute the subtract-off gamma1 for D1
# Compute the pihat subtract-off by projecting the influences onto the scores.
# Pseudoinverse is important if some combo of J*Delta is always/never observed.
gamma1 <- scores %*% ginv(t(scores)%*%scores) %*% (t(scores) %*% D1)

# This is the correct variance of the baseline cumhaz at the time of interest.
# This accounts for the estimated weights.
D3pihat <- (D1-gamma1) - ( (D2-gamma2) %*% EdLambda )
varBaselineCumhaz <- t(D3pihat) %*% D3pihat
#print( varBaselineCumhaz )

# Now get the covariance of the baseline cumhaz and beta
covLambdaBeta <- t(D3pihat) %*% D2pihat  


# Create variance matrix of all 768 possible covariate-specific survivals that are:
# Z = 2smoke*2drink*2bmi*6agestrsex*4selquartiles*4viteadjquartiles
# and print out the survival for each of those 768 covariates.
#
# expand.grid makes a single replicate design matrix for all these, making sure to label
# them as factors with as many levels as each variable.  Since the left column varies the
# most, I reverse the matrix so that the column varying the most is on the right (this
# must be the vitamin e column since that's the last variable in the coxmod).  Once it's
# reversed back with rev, the covars are in the same order as coxmod.  It returns the
# variables named as Var1,Var2,etc.  Then model.matrix expands out the factors into a
# model matrix using the default contrasts (contr.treatment).  The model formula gets rid
# of the intercept (there is none with the coxmod) and the variable names appear reversed
# because I had to reverse earlier.  Finally, the first column is redundant, so I get rid
# of it.

#
# Some partial sample output:
# From expand.grid:
#     Var1 Var2 Var3 Var4 Var5 Var6
# 1      1    5    9   15   17   19
# 2      2    5    9   15   17   19
# 3      3    5    9   15   17   19
# 4      4    5    9   15   17   19
# 5      1    6    9   15   17   19
# 6      2    6    9   15   17   19
# 7      3    6    9   15   17   19
# 8      4    6    9   15   17   19
# This gives me all possible covariates (a 1 replicate design matrix)
#
# From model.matrix:
#     Smoke- Smoke+ drink+   BMI+ femmid femold malyng malmid malold  sel2  sel3  sel4 vite2 vite3 vite4
#     Var619 Var620 Var518 Var416 Var310 Var311 Var312 Var313 Var314 Var26 Var27 Var28 Var12 Var13 Var14
# 1        1      0      0      0      0      0      0      0      0     0     0     0     0     0     0
# 2        1      0      0      0      0      0      0      0      0     0     0     0     1     0     0
# 3        1      0      0      0      0      0      0      0      0     0     0     0     0     1     0
# 4        1      0      0      0      0      0      0      0      0     0     0     0     0     0     1
# 5        1      0      0      0      0      0      0      0      0     1     0     0     0     0     0
# 6        1      0      0      0      0      0      0      0      0     1     0     0     1     0     0
# 7        1      0      0      0      0      0      0      0      0     1     0     0     0     1     0
# 8        1      0      0      0      0      0      0      0      0     1     0     0     0     0     1
# The first 2 columns are for smoking which is binary, this is redundant, so I get rid of
# 1st column.
# This is the prototype code that I am emulating:
# Z <- model.matrix(~ -1+Var6+Var5+Var4+Var3+Var2+Var1,
#                   data.frame(rev(expand.grid(factor(1:4),factor(5:8),factor(9:14),factor(15:16),
#                                              factor(17:18),factor(19:20)))))

# Emulate the above code as follows: Get each covariate in a list.  This is in the
# "variables" attribute of the terms object.  This attribute is a call object, so
# evaluating it leaves us with each covariate in a list.  The first member of the list is
# the reponse, so strip it off.
covarlist <- eval(attr(coxmod$terms,"variables"),envir=data)[-1]

# Generates the model frame (model matrix but no aliasing or dummying) I get the number of
# levels for each factor, make a list of 1:nlevels for each factor, generate the model
# frame, reverse it, convert to each to factor (reversing is needed to get the factor that
# changes the most to be on the rightmost column, that's just the way I like it).  Note
# that the S function fac.design() doesn't exist in R, but you can approximate it with:
# expand.grid(lapply(list(2,3,2),seq,from=1)), makes 1 replicate of a full 2x3x2
# factorial.
temp <- as.data.frame(lapply(rev(expand.grid(rev(lapply(lapply(covarlist,nlevels),
                                                        seq,from=1)))),as.factor))

# Generate the model matrix Generate the formula as a string of the sum of each of the
# variable names, without the intercept For example if I have 3 variables: "~
# -1+Var3+Var2+Var1". Then evaluate this with as.formula()
Z <- model.matrix(as.formula(paste(c("~ -1",names(temp)),collapse="+")),data=temp)

# Now I've emulated the code to make the full model matrix.  The first column is
# unnecessary
Z <- Z[,-1]

ncovars <- dim(Z)[1]
#print(covariatesurvs <- exp(-cumhaz*exp(Z %*% beta)))
covariatesurvs <- exp(-cumhaz*exp(Z %*% beta))
ncolDpihat = 1+length(beta)
Dpihat <- matrix(c(D3pihat,D2pihat),ncol=ncolDpihat,byrow=FALSE)
V <- t(Dpihat) %*% Dpihat
G <- matrix(exp(-cumhaz*exp(Z %*% beta))*exp(Z %*% beta),nrow=ncovars,ncol=ncolDpihat) *
     matrix(c(rep(1,ncovars),cumhaz*Z),nrow=ncovars,ncol=ncolDpihat)
#print( varmat <- G %*% V %*% t(G) )
varmat <- G %*% V %*% t(G)



##
# Now create adjusted survivals and the population-marginal crude survival for the entire
# cohort.  Also get their variances and covariances.  Get adjusted survivals for each
# level of the exposures.
##

# First get the matrices for all the exposures, and get the number of levels
exposurelist <- eval(attr(terms(as.formula(paste(c(outcome,"~",exposures),collapse=""))),
                          "variables"),envir=data)[-1]
exposurelevels <- unlist(lapply(exposurelist,nlevels))
nexposures <- prod(exposurelevels)


# First get the probabilities of the confounders, which are
# smoke+drink+bmimedian+agestr*sex Note that bmi,smoke and drink are missing on a few
# people, so assume they are MCAR missing and and just divide by the total of people who
# have all observed.  These are used to compute the adjusted survival.  I need to compute
# the number of people in all cells of the confounder combinations:
# Ex: if confounders are
# "smoke+drink+interaction(agestr,sex)" then the command I want is:
# padj <- table(interaction(interaction(agestr,sex),drink,smoke))
if (confounders=="") {
  # There are no confounders, so nothing to adjust for
  padj <- 1
}
else {
  # Compute the joint confounder distribution
  padj <- paste(c("table(interaction(",
                  paste(c(rev(unlist(strsplit(confounders,"\\+")))),collapse=",")
                  ,"))"),collapse="")
  padj <- eval(parse(text=padj),envir=data)
  padj <- padj/sum(padj)
}

# Next get the probabilities of the entire covariate distribution, by summing the weights
# for each covariate combination amongst those with the entirely observed covariate
# vector.  Any combination that is NA means zero were observed, so I set them manually to
# zero.  This is used to compute the crude survival Need to get weighted counts in each
# cell for all covariate combinations:
# Ex: if confounders are "smoke+drink+bmimedian+interaction(agestr,sex)" and
# if exposures are "selquartiles+viteadjquartiles", then the command I want is:
# pcrude <- tapply(1/p.i.h.a.t.[o.b.s.e.r.v.e.d.==TRUE],
#                  interaction(viteadjquartiles,selquartiles,agestr,sex,bmimedian,drink,smoke)[o.b.s.e.r.v.e.d.==TRUE],
#                  sum)
if (confounders=="") {
  # There are no confounders, so just do the exposures
  pcrude<-paste(c("tapply(1/p.i.h.a.t.[o.b.s.e.r.v.e.d.==TRUE],interaction(" ,
                  paste(c(rev(unlist(strsplit(exposures,"\\+")))),collapse=",")
                  , ")[o.b.s.e.r.v.e.d.==TRUE],sum)"),collapse="")
}
else {
  # Get joint distribution of exposures and confounders
  pcrude<-paste(c("tapply(1/p.i.h.a.t.[o.b.s.e.r.v.e.d.==TRUE],interaction(" ,
               paste(c(rev(unlist(strsplit(paste(c(confounders,exposures),collapse="+"),"\\+")))),collapse=",")
                  , ")[o.b.s.e.r.v.e.d.==TRUE],sum)"),collapse="")
}
pcrude <- eval(parse(text=pcrude),envir=data)
pcrude[is.na(pcrude)] <- 0
normalization <- sum(pcrude)
pcrude <- pcrude/sum(pcrude)

# The contrast matrix for the 3 adjusted survivals are just 48 matrices each of whose
# elements are zero or a single element of padj, linked in columns (so is 3 x
# ncovars). Each of the 16 rows are one of the 16 selenium*vitamine combinations, the 3
# rows are for Q1Q1, QmQm and Q4Q4 respectively.  Within each block we need the
# probability of the confounder space constant, so we need one 1 such block for each level
# of the confounders.
#C <- NULL
#for (i in 1:48) { C <- cbind(C,diag(padj[i],nrow=16,ncol=16)) }
C <- kronecker(t(padj),diag(nexposures))

# The contrast matrix for the crude survival is just the covariate combination specific to
# each covariate survival.  So the contrast matrix is (nexposures+1) x ncovars
C <- rbind(C,pcrude)

# Adjusted Survivals and Crude Survival (crude survival is the last element of vector)
allsurvs <- C %*% covariatesurvs
# Their variances and covariances
varallsurvs <- C %*% varmat %*% t(C)
stderrallsurvs <- sqrt(diag(varallsurvs))
exposurenames <- gsub("[+]",",",exposures)
# Get CI via complementary-loglog transformation
theta <- exp(-1.96*stderrallsurvs / log(allsurvs^allsurvs))

# If timeofinterest is Inf, set to end of followup (this was set to last riskset time)
if (is.infinite(usertimeofinterest)) {
  timeofinterest <- max(lifetime) # Print out end of followup for user
}

# Output the table, either on the cumulative incidence scale or survival scale
if ( cuminc ) {
  survtable <- cbind(1-allsurvs,stderrallsurvs,1-allsurvs^(1/theta),1-allsurvs^theta)
  colnames(survtable) <- c("Cumulative Incidence","StdErr","95% CI Left","95% CI Right")  
  cat("\n")
  cat(c("Standardized Cumulative Incidence for",exposurenames,"by time",timeofinterest,
        "\n"))
}  
else {
  survtable <- cbind(allsurvs,stderrallsurvs,allsurvs^theta,allsurvs^(1/theta))
  colnames(survtable) <- c("Survival","StdErr","95% CI Left","95% CI Right")
  cat("\n")
  cat(c("Standardized Survival for",exposurenames,"by time",timeofinterest,"\n"))
}

# Create the rownames: Get levels for each exposure, make a full replicate of that
# (reverse it as I prefer to have the fastest changing level on the right), convert each
# cell to a character, then paste the columns together with a comma in between, put in the
# exposure variable names, add "Crude" as the last row.
# Ex: rownames(survtable) <- c("Sel1 VitE1","Selmiddle VitEmiddle","Sel4 VitE4","Crude")
# Example that works if you only have 2 exposures:
# rownames(survtable) <- paste(exposurenames,as.vector(sapply(levels(temp),paste,levels(temp),sep=",")))
exposurereplicate <- apply(apply(rev(expand.grid(lapply(exposurelist,levels))),
                                 c(1,2),as.character),
                           1,paste,collapse=",")
# If you want the names of the exposure variables repeated in each row, do this:
#rownames <- paste(exposurenames,exposurereplicate,sep=": ")
rownames <- exposurereplicate
rownames(survtable) <- c(rownames,"Crude")

# Output Survival table
print(survtable)


##
# Use the exposureofinterest to compute risk differences and PAR towards
# If none is given by the user, set to the first level of the exposures
# Use partial matching to find the unique exposure from the user's input, if none exists 
# or the match isn't unique, then return an error.
##
if (is.na(exposureindex <- pmatch(exposureofinterest,exposurereplicate))) {
  if (exposureofinterest=="") {
    exposureindex <- 1
  }
  else {
    cat("\n The ",length(exposurereplicate),
        "possible exposures of interest that you can specify are:\n")
    print(exposurereplicate)
    stop("The entered exposure of interest doesn't match any in the above list")
  }
}
# Use the name of what this program matched to the user's exposureofinterest, so that user
# knows what his argument got matched to, since this is what is used by this program.
exposureofinterest <- exposurereplicate[exposureindex]

# Compute the risk differences towards the exposureofinterest between standardized
# survivals
riskdiffs <- allsurvs[exposureindex]-allsurvs[-exposureindex]

# The contrast matrix for the variances is diagonal with -1 in diagonals, splicing inside
# a column where the exposure of interest is, that is a column of ones.  Note that this
# includes the crude survival in a contrast, so the matrix is nexposures X nexposures+1.
# I make the diagonal matrix, cut it, splice in the 1s, and glue it back together.
contrasts <- -diag(nexposures)
contrasts <- cbind(contrasts[,1:(exposureindex-1)],1,contrasts[,-(1:(exposureindex-1))])

# Compute standard errors and CI and output in table
stderrriskdiffs <- sqrt(diag(contrasts %*% varallsurvs %*% t(contrasts)))
CIriskdiffs <- cbind(riskdiffs-1.96*stderrriskdiffs,riskdiffs+1.96*stderrriskdiffs)
riskdifftable <- t(rbind(riskdiffs,stderrriskdiffs,t(CIriskdiffs)))
colnames(riskdifftable) <- c("Risk Difference","StdErr","95% CI Left","95% CI Right")
rownames(riskdifftable) <- paste(exposureofinterest,"-",rownames(survtable)[-exposureindex])
cat("\n")
cat(c("Standardized Risk Differences vs.",exposurenames,"=",exposureofinterest,
      "by time",timeofinterest,"\n"))
print(riskdifftable)

##
# Compute PAR = 1 - P(D|~E)/P(D) 
# The proportion of disease eradicated if everyone has the disease rate of those in the
# 4th quartiles of selenum and cholesterol-adjusted Vitamin E.
# The numerator is this adjusted disease rate, the denominator is the crude disease rate
#
# Also compute Risk Difference of quartile 4 to the crude
##
adjsurv <- allsurvs[exposureindex]
crudesurv <- allsurvs[nexposures+1]
varadjandcrude <- varallsurvs[c(exposureindex,nexposures+1),c(exposureindex,nexposures+1)]
PAR <- 1 - (1-adjsurv)/(1-crudesurv)

# Compute the variance of PAR and it's CI based on the log(1-PAR) transform (I call it
# logPAR)
myD <- matrix( c( 1/(1-crudesurv), (adjsurv-1)/(1-crudesurv)^2 ), nrow=2,ncol=1)
myvarPAR <- t(myD) %*% varadjandcrude %*% myD
myCIPAR <- c(PAR-1.96*sqrt(myvarPAR),PAR+1.96*sqrt(myvarPAR))
D <- (C[exposureindex,]-C[nexposures+1,])*(1-crudesurv)^-1 +
    ((C[exposureindex,]-C[nexposures+1,])%*%covariatesurvs)*(1-crudesurv)^-2*C[nexposures+1,]
varPAR <- D %*% varmat %*% t(t(D))
CIPAR <- c(PAR-1.96*sqrt(varPAR),PAR+1.96*sqrt(varPAR))
varlogPAR <- (-1/(1-PAR))^2 * varPAR
CIlogPAR <- 1-exp(c(log(1-PAR)+1.96*sqrt(varlogPAR),log(1-PAR)-1.96*sqrt(varlogPAR)))

PARtable <- t(as.matrix(c(PAR,sqrt(varPAR),CIlogPAR)))
colnames(PARtable) <- c("Estimate","StdErr","95% PAR CI Left","95% PAR CI Right")
rownames(PARtable) <- "PAR"
cat("\n")
cat("PAR if everyone had",exposurenames,"=",exposureofinterest,"\n")
print(PARtable)
cat("\n")

##
# If we want, plot the adjusted survivals
##
if ( plot==TRUE ) {

  ##
  # Compute the entire Adjusted Survival Curve for each of the exposures + 1 for crude
  # (Not that I use the crude)
  ##
  ntimes <- length(time)
  if ( cuminc )
    adjsurvcurv <- matrix(0,nrow=nexposures+1,ncol=ntimes)
  else
    adjsurvcurv <- matrix(1,nrow=nexposures+1,ncol=ntimes)
  cumhazs <- cumsum(hazard)
  for (timeindex in 2:ntimes) {
    cumhaz <- cumhazs[timeindex-1]
    covariatesurvs <- exp(-cumhaz*exp(Z %*% beta))
    if ( cuminc )
      adjsurvcurv[,timeindex] <- 1 - C %*% covariatesurvs 
    else
      adjsurvcurv[,timeindex] <- C %*% covariatesurvs
  }
  
  # Plot cumulative incidence curves
  # First setup the plot
  # The smallest y is the smallest survival or the lowest lower end of the CI, whichever
  # is smaller.  Vice versa for cumulative incidence.

  if ( plotfilename != "" )
    pdf(plotfilename)

  if ( cuminc )
    plot(time,y=NULL,type="n",ylim=c(0,max(max(adjsurvcurv),max(survtable[,4]))),
         xlim=c(0,max(time)),xlab="Time",ylab="Standardized Cumulative Incidence",...)
  else
    plot(time,y=NULL,type="n",ylim=c(min(min(adjsurvcurv),min(survtable[,3])),1),
         xlim=c(0,max(time)),xlab="Time",ylab="Standardized Survival",...)
  

  endofstudy <- length(time)
  for ( i in 1:nexposures ) {

    # Plot the curve
    lines(time,adjsurvcurv[i,],type="s")

    # Put in error bars at the time of interest
    # step is the size of the horizontal bar on top/bottom, length set to 2% of the
    # maximum x-axis value
    lines(rep(time[risksetofinterest],2),survtable[i,3:4])
    step <- max(time)*0.01
    lines(c(time[risksetofinterest]-step,time[risksetofinterest]+step),
          rep(survtable[i,3],2))
    lines(c(time[risksetofinterest]-step,time[risksetofinterest]+step),
          rep(survtable[i,4],2))

    # Notate each curve
    text(time[endofstudy]-0*step,adjsurvcurv[i,endofstudy],pos=1,exposurereplicate[i])

  }

  if ( plotfilename != "" )
    dev.off()
}




# If no plotting for stairs available, use the kronecker product to create all the x,y
# pairs for a survival curve
#plot(c(0,kronecker(t(time[-1]),c(1,1))),c(kronecker(t(adjsurvcurv[1,-endofstudy]),c(1,1)),adjsurvcurv[1,endofstudy]),type="l")

# Return this
output <- list()
output$coxmod <- coxmod
output$samplingmod <- samplingmod
output$survtable <- survtable
output$riskdifftable <- riskdifftable
output$PARtable <- PARtable
if ( plot==TRUE ) {
  # Output the time, the adjusted survivals, and the crude survival.
  # In the matrix, make sure the row names can be legitimate R variables.
  output$plotdata <- cbind(time,t(adjsurvcurv))
  colnames(output$plotdata) <- make.names(c("time",rownames(survtable)),unique=TRUE)
}

return(output)
}
