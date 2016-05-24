nested.km <-
  ##
  # Compute Kaplan-Meier non-parametric survival within groups defined by strata
  # using pihat weights for studies nested within cohorts.
  # By: Hormuzd Katki 4/10/09
  #
  # DO NOT USE ANY STRATA() STATEMENTS (which is illegal for survfit())
  # No left truncation, no late entry allowed -- everyone must enter at a chosen time zero
  # Reason why covariate is missing cannot depend on time (no incidence-density sampling).
  # Does not account for competing risks.
  # Requires a data frame, name no variables in it 'o.b.s.e.r.v.e.d.' or 'p.i.h.a.t.'
  ##
  function(survfitformula, samplingmod, data, outputsamplingmod=FALSE,outputriskdiff=FALSE,
           exposureofinterest = "", timeofinterest = Inf, glmlink = binomial(link="logit"),
           glmcontrol = glm.control(epsilon=1e-10, maxit=10, trace=FALSE),
           missvarwarn = TRUE,...)
{

# Need generalized inverses from MASS
library(MASS)

# Check if data.frame is entered
if ( is.null(data) | !is.data.frame(data) )
  stop("data is a required argument and must be a data.frame")

# Check that there are no variables in the dataset named "o.b.s.e.r.v.e.d."
if ( any("o.b.s.e.r.v.e.d." == colnames(data)) )
  stop("You are not allowed to name any variable in samplingmod as 'o.b.s.e.r.v.e.d.'.")

# Check that there are no variables named "p.i.h.a.t."
if ( any("p.i.h.a.t." == colnames(data)) )
  stop("You are not allowed to name any variable in the data as 'p.i.h.a.t.'.")

# Make sure that time of interest is a number.  
if (!is.numeric(timeofinterest) | timeofinterest<0)
  stop("Time of interest must be between 0 and end of follow-up (use Inf to mean end of follow-up)")

# Make sure strata have no '*', or strata(),cluster(),offset() statements
illegal <- c("*","strata(","cluster(","offset(")
if ( any(unlist(lapply(illegal,grep,x=survfitformula,fixed = TRUE))) )
  stop(paste("The survfitformula contains some illegal commands:", paste(illegal,collapse=", ")))

# Extract the Surv() object out of the model formula, keep only the observed
# First column is the lifetime, second column the event indicator
# Make sure everyone enters cohort simultaneously
outcomelist <- eval(attr(terms(as.formula(survfitformula)),"variables"),envir=data)[1]
Survdata <- outcomelist[[1]]
if ( dim(Survdata)[2] > 2 )
  stop("Left-truncation (late entry) not supported: all subjects must enter cohort simultaneously")

# Get the strata, which is a list of the covariates in the formula
covarlist <- eval(attr(terms(as.formula(survfitformula)),"variables"),envir=data)[-1]

# Generate vector of which observations have a fully-observed covariate vector
# This is 1 if the stratum is observed for that subject, 0 otherwise
# Then make the formula for the sampling model
observed <- as.numeric(complete.cases(covarlist))
# Put observed into the dataframe as o.b.s.e.r.v.e.d.
data$o.b.s.e.r.v.e.d. <- observed
if (samplingmod != "")
  samplingformula <- paste("o.b.s.e.r.v.e.d. ~",samplingmod)
else # If nothing to stratify on, then just regress on the intercept
  samplingformula <- paste("o.b.s.e.r.v.e.d. ~ 1")

# Make sure that strata are factors
if ( any(sapply(covarlist,inherits,"factor")==FALSE) )
  stop("All strata must inherit from class factor")

# If you are missing data in the samplingmod, get rid of it, but issue a warning
temp <- cbind(outcomelist,
              as.data.frame(eval(attr(terms(as.formula(samplingformula)),"variables"),
                                 envir=data)))
missing <- apply(is.na(temp),1,any)
if ( any(missing) ) {
  if ( missvarwarn) {
    warning(paste("You had",sum(missing),
                  "obervations with missing values in the samplingmod",
                  "or survival outcome.\n",
                  "These observations were removed.\n",
                  "There shouldn't be any missing data in the samplingmod or survival",
                  "outcome, so be careful.\n"
                  ))
  }
  data <- data[!missing,]

  # Regenerate the covarlist, survival outcome, and who's observed:
  # Generate vector of which observations have a fully-observed covariate vector
  # This is 1 if all exposures and confounders are observed, 0 otherwise
  # Note that I don't need to regenerate o.b.s.e.r.v.e.d. since that's taken care of in
  # covarlist.  Also regenerate the survival outcome.
  covarlist <- eval(attr(terms(as.formula(survfitformula)),"variables"),envir=data)[-1]
  observed <- as.numeric(complete.cases(covarlist))
  Survdata <- eval(attr(terms(as.formula(survfitformula)),"variables"),envir=data)[[1]]

}

# Cohort size
n <- dim(data)[1]


# I must use treatment and polynomial contrasts, so set these temporarily, then revert
# back to the user's settings.  I got this code from the Package haplo.stats()
contrasts.tmp <- c(factor = "contr.treatment", ordered = "contr.poly")
contrasts.old <- options()$contrasts
options(contrasts = contrasts.tmp)
on.exit(options(contrasts = contrasts.old))

# First fit sampling model; remember to keep the model matrix by setting x=T.
# For ease, use treatment and polynomial contrasts
# Make sure that no covariate combination has zero chance of being observed
suppressWarnings(
                 samplingmod <- glm(as.formula(samplingformula), x=TRUE, family=glmlink, 
                                    na.action=na.fail, data=data, control=glmcontrol,...)
                 )

# Get the estimated pihats to use for weighting.
pihat <- data$p.i.h.a.t. <- predict(samplingmod,type='response')

# If any pihats are exactly zero, stop right here -- cannot make estimates!
if ( any(pihat==0) )
  stop(paste(
             "According to the sampling model, some subjects have zero probability of",
             "being observed. \n",
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

# Fit the Non-parametric model with weights, to all strata.  The estimates will be
# correct, but their standard errors will be totally wrong; all of the lines after this
# will produce the correct standard error
survmod <- survfit(as.formula(survfitformula), type="fl",weights=1/pihat,
                   na.action=na.omit,data=data,...)

# Only keep those who are not missing lifetime or event indicator
keep <- ( !is.na(Survdata[,1]) & !is.na(Survdata[,2]) )
Survdata <- Survdata[keep, ]
data <- data[keep,]
covarlist <- as.data.frame(covarlist)
covarlist <- covarlist[keep,]

# Figure out which is the stratum (exposure) of interest
# If none is given by the user, set to the first level of the exposures
# Use partial matching to find the unique exposure from the user's input, if none exists 
# or the match isn't unique, then return an error.
# First make sure covarlist is really a list (it won't be if there's only one exposure variable)
if ( !is.list(covarlist) )
  covarlist <- list(covarlist)
allexposures <- apply(apply(rev(expand.grid(lapply(covarlist,levels))), c(1,2),
                            as.character), 1, paste, collapse=",")
if (is.na(exposureindex <- pmatch(exposureofinterest,allexposures))) {
  if (exposureofinterest=="") {
    exposureindex <- 1
  }
  else {
    cat("\n The ",length(allexposures),"possible exposures of interest that you can specify are:\n")
    print(allexposures)
    stop("The entered exposure of interest doesn't match any in the above list")
  }
}

# Use the name of what this program matched to the user's exposureofinterest, so that user
# knows what his argument got matched to, since this is what is used by this program.
# Using survmod$strata includes the names of the variables on the data frame that were
# used.
exposureofinterest <- names(survmod$strata)[exposureindex]


# This line breaks the dataframe data into data.frames (by applying the function
# data.frame) according to all combinations of levels of the factors in covarlist.  These
# combinations are the strata.  The factors are reversed so that the order of the
# combinations is equivalent to the order of the strata in survmod.  Only want to keep
# those with observed data, but this is automatically done because NA is not a level in
# the exposure factors, so any missing exposures are dropped -- this could be a problem if
# the levels of the exposure factors includes NA.
strata <- by(data,rev(covarlist),data.frame)
# Keep only non-empty strata
nullstrata <- sapply(strata,is.null)
strata <- strata[which(nullstrata==FALSE)]

# Do the same for Survdata
Survstrata <- by(Survdata,rev(covarlist),data.frame)
# Keep only non-empty strata
Survstrata <- Survstrata[which(nullstrata==FALSE)]

# Get the interaction of all the stratum-defining factors, so we can denote the stratum
# each row belongs to.  Reverse the covariate list to order the strata in the same order
# as survfit().  Unclass to create integer numeric codes for each level, these go from 1
# to #strata.  Note that this forces to NA any row with a missing covariate, so this
# indicates on those with observed=1, which is what I want.  But if NA is defined as a
# level of the exposure covariates, they would mistakenly count as another stratum, and
# this might screw up.  Note that I must say drop=FALSE, even though I want drop=T, because
# that changes the ordering of the factor levels!
stratindic <- unclass(interaction(rev(covarlist),drop=FALSE))

# Drop any empty strata (since drop=T can't be used).  Find the strata that are not in the
# set of nullstrata, keep them along with any subjects with missing exposures.  We want
# the kept strata to have consecutive numbering, which won't be the case if we remove the
# null strata, so redo the numbering by converting back to factor and then unclassing to
# create integer numeric codes for each level from 1 to #strata.
stratindic <- unclass(factor(stratindic[is.element(stratindic,
                                                   c(NA,NaN,which(nullstrata==FALSE)))]))

##
# Loop through each time in each stratum
##
nstrata <- length(strata)
risksetofinterest <- rep(0,nstrata)
keepsurvtablerow <- matrix(0,nrow=nstrata,ncol=5)
pihatinfluences <- matrix(0,nrow=n, ncol=nstrata)
for ( s in 1:nstrata ) {

  # Get data for the current stratum, by def'n, this gets rid of anyone missing exposures
  stratum <- strata[[s]]
  Survstratum <- Survstrata[[s]]
  hpstrat <- survmod[s]
  thisstratum <- (stratindic==s)
  thisstratum[is.na(thisstratum)] <- FALSE # Anyone with missing exposures is in no strata
  
  # Only keep the event times (not censored times) in that stratum.
  # Extract out the time, S0, and weighted-dNbar, and remember to add in time=0,
  # S0(time=0), dNbar(time=0)=0.
  failtimes <- hpstrat$n.event!=0
  time <- c(0,hpstrat$time[failtimes])
  S0 <- c(sum(1/stratum$p.i.h.a.t.),hpstrat$n.risk[failtimes])
  dNbar <- c(0,hpstrat$n.event[failtimes])

  # Check if timeofinterest is a legal time, if it's Inf (or too big), set to max(time)
  # Get the risket number of the time of interest
  # Don't change the user's timeofinterest (unless it's Inf), remember that for the output
  ourtimeofinterest <- timeofinterest
  if (is.infinite(ourtimeofinterest) | timeofinterest>max(time)) {
    ourtimeofinterest=max(time)
    timeofinterest <- max(hpstrat$time) # Print out end of followup for user
  }
  risksetofinterest[s] <- sum(time<=ourtimeofinterest) 

  # Compute riskset (the index of the last failure time at risk) for each subject.
  nsubjects <- dim(stratum)[1]
  riskset <- matrix(0,nrow=nsubjects,ncol=1)
  
  # This below line works in R, but SPLUS doesn't have the which.max function -- which
  # finds the the index of the maximum of a vector.  If there are ties, it takes the first
  # index.  I want the last index which is why I have to use rev to reverse the vector.
  # Use the below 2 lines in SPLUS, workaround for not having the
  # which.max function
  #   for (i in 1:nsubjects) {
  #     eligible <- rev(time <= Survstratum[i,1])
  #     riskset[i] <- (length(time)+1)-which(eligible==max(eligible))[1]
  #   }

  for (i in 1:dim(stratum)[1]) {
    riskset[i] <- (length(time)+1)-which.max(rev(time <= Survstratum[i,1]))
  }
  
  
  # Compute Influence Function for all observed subjects at each failure time by looping
  # through each timeindex r
  correction <- cumsum(S0^-2 * dNbar)
  ntimes <- length(time)-1
  influence <- matrix(0,nrow=nsubjects,ncol=ntimes)
  for (r in 2:(ntimes+1)) {
    influence[,r-1] <-(1/stratum$p.i.h.a.t.)*
                      ((1/S0[riskset])*(Survstratum[,2]==1)*(riskset<=r)
                                       -  correction[pmin(riskset,r)] )
  }
  
  # Give all those with missing covariates a zero influence
  D1 <- matrix(0,nrow=n,ncol=ntimes)
  D1[thisstratum,] <- influence

  # Compute the pihat subtract-off by projecting the influences onto the scores from
  # logistic model.  Pseudoinverse is important if some combo of the samplingmod are
  # always/never observed.
  gamma1 <- scores %*% ginv(t(scores)%*%scores) %*% (t(scores) %*% D1)

  # Keep track of the pihat influence functions for each stratum.  We'll need this to
  # compute variances of risk differences between strata.  (Ignore the zeroth riskset, so
  # subtract 1)
  pihatinfluences[,s] <- (D1-gamma1)[,risksetofinterest[s]-1]

  # This is the correct variance of the cumhaz in the each stratum at all times.
  varcumhaz <- diag(t(D1-gamma1)%*%(D1-gamma1))

  # Delta-Method yields the variance of the survivals
  surv <- hpstrat$surv[hpstrat$n.event!=0]
  varsurv <- surv^2 * varcumhaz

  # CI on cloglog scale
  stderrcloglog <- sqrt(log(surv)^-2 * varcumhaz)
  criticalvalue <- qnorm((1-survmod$conf.int)/2)
  CI <- cbind( surv^(exp(-criticalvalue*stderrcloglog)) ,
               surv^(exp(+criticalvalue*stderrcloglog)) )
    
  # Print out the table (suppressed)
  survtable <- cbind(time[-1],hpstrat$surv[hpstrat$n.event!=0],sqrt(varsurv),CI)
  keepsurvtablerow[s,] <- survtable[risksetofinterest[s]-1,]
  colnames(survtable) <- c("Time","Survival","Std.Err","Left CI","Right CI")

  # Put inside the survfit object
  # First figure out the indices in the $surv, $std.err, $upper, $lower where this
  # stratum is
  beginingindex <- cumsum(c(1,survmod$strata))
  left <- beginingindex[s]
  indices <- left:(beginingindex[s+1]-1)
  # Put in the survtable stderrs and CI at the failure times
  failureindices <- (left-1)+which(hpstrat$n.event>0)
  # not stderr of survival survtable[,3], but summary(survmod) returns std.err of survival
  survmod$std.err[failureindices] <- sqrt(varcumhaz) 
  survmod$lower[failureindices]   <- survtable[,4]
  survmod$upper[failureindices]   <- survtable[,5]
  # Now update the stderrs and CI at the censoring times
  censorindices <- (left-1)+which(hpstrat$n.event==0)
  # The censoring times get the values at the previous time, except # if survival=1.
  # Fix this up later in this function
  copyindices <- censorindices-1 
  copyindices[copyindices==(left-1)] <- NA 

  # Assign the censoring times the values at the previous time (but we really want this to
  # be the previous failure time) Since the previous time might not be a failure (in case
  # of multiple censorings in a row) we have to loop through each time sequentially to
  # that the second censoring gets the value from the first censoring which got it from
  # the previous failure.
  if ( length(censorindices)>0 ) 
    for (i in 1:length(censorindices)) {
      currentcensorindex <- censorindices[i]
      currentcopyindex <- copyindices[i]
      # If any censored before the first failure, their survival is 1, stderr=0, CI=(1,1)
      if ( is.na(currentcopyindex) ) {
        survmod$std.err[left] <- 0
        survmod$lower[left] <- 1
        survmod$upper[left] <- 1
      }
      else { # Assign the values at the current censoring time to that of the previous time
        survmod$std.err[currentcensorindex] <- survmod$std.err[currentcopyindex]
        survmod$lower[currentcensorindex]   <- survmod$lower[currentcopyindex]  
        survmod$upper[currentcensorindex]   <- survmod$upper[currentcopyindex]
      }
    }                    
  
} # End looping through all strata and event-times


##
# Now print risk differences vs. exposure of interest and CIs at the time of interest
##

# Compute the risk differences against the exposureofinterest
riskdiffs <- keepsurvtablerow[exposureindex,2] - keepsurvtablerow[-exposureindex,2]

# Compute the variance-covariance matrix of the survivals
cumhazvarmat <- t(pihatinfluences) %*% pihatinfluences
survs <- matrix( rep(keepsurvtablerow[,2],nstrata) , nrow=nstrata,byrow=TRUE)
survvarmat <- survs^2 * cumhazvarmat

# The contrast matrix for the variances is diagonal with -1 in diagonals, splicing inside
# a column where the exposure of interest is, that is a column of ones.  The matrix is
# nstrata X nstrata.  I make the diagonal matrix, put in the 1s for the column
# representing the exposureofinterest, and get rid of the row that represents just the
# exposureofinterest.
contrasts <- -diag(nstrata)
contrasts[,exposureindex] <- 1
contrasts <- contrasts[-exposureindex,]
if ( is.vector(contrasts) )
  contrasts <- t(as.matrix(contrasts))

# Compute standard errors and CI and output in table
stderrriskdiffs <- sqrt(diag(contrasts %*% survvarmat %*% t(contrasts)))
CIriskdiffs <- cbind(riskdiffs-1.96*stderrriskdiffs,riskdiffs+1.96*stderrriskdiffs)
riskdifftable <- t(rbind(riskdiffs,stderrriskdiffs,t(CIriskdiffs)))
colnames(riskdifftable) <- c("Risk Difference","StdErr","95% CI Left","95% CI Right")
rownames(riskdifftable) <- paste(allexposures[exposureindex],"-",
                                 allexposures[-exposureindex])

# Get the variable(s) name(s) that were used to make the strata
exposurenames <- unlist(strsplit(names(survmod$strata)[1],split="=",fixed=TRUE))[1]
cat("\n")
cat(c("Risk Differences vs.",exposureofinterest,"by time",timeofinterest,"\n"))
print(riskdifftable)
cat("\n")
    
# Return just survfit model, sampling models, or risk difference table
output <- survmod
if ( outputsamplingmod | outputriskdiff ) {
  output <- list(survmod=survmod)
  if ( outputsamplingmod )
    output$samplingmod <- samplingmod
  if ( outputriskdiff )
    output$riskdiff <- riskdifftable
}
return(output)

} # End function
