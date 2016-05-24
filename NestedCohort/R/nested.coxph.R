nested.coxph <-
  ##
  # Compute Cox model hazard ratios from studies nested within cohorts
  # By: Hormuzd Katki 4/10/09
  ##
  function(coxformula, samplingmod, data, outputsamplingmod = FALSE,
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
  
# Generate vector of which observations have a fully-observed covariate vector
# This is 1 if all exposures and confounders are observed, 0 otherwise
# Then make the formula for the sampling model
covarlist <- eval(attr(terms(as.formula(coxformula)),"variables"),envir=data)[-1]
outcome <- eval(attr(terms(as.formula(coxformula)),"variables"),envir=data)[1]
observed <- as.numeric(complete.cases(covarlist))

# Make sure cox model formula has no cluster() statements
illegal <- c("cluster(")
if ( any(unlist(lapply(illegal,grep,x=coxformula,fixed=TRUE))) )
  stop(paste("exposures or confounders contain some illegal commands:",
             paste(illegal,collapse=", "),"\n"))

# Put observed into the dataframe and create the sampling formula
data$o.b.s.e.r.v.e.d. <- observed
samplingformula <- paste("o.b.s.e.r.v.e.d. ~",samplingmod)

# Check that there are no variables in the dataset named "o.b.s.e.r.v.e.d."
if ( any("o.b.s.e.r.v.e.d." == all.names(as.formula(samplingformula))[-2]) )
  stop("You are not allowed to name any variable in samplingmod as 'o.b.s.e.r.v.e.d.'.")

# Check that there are no variables named "p.i.h.a.t."
if ( any("p.i.h.a.t." == colnames(data)) )
  stop("You are not allowed to name any variable in the data as 'p.i.h.a.t.'.")

# If you are missing data in the samplingmod or outcome, get rid of it,
# but issue a warning to the user
temp <- cbind(outcome,
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
  # This is 1 if all exposures and confounders are observed, 0 otherwise
  # Then make the formula for the sampling model
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

# First fit sampling model; remember to keep the model matrix by setting x=TRUE.  Make sure
# that no covariate combination has zero chance of being observed.  Suppress warnings that
# model hasn't converged -- this is usually because some fitted probabilities are close to
# 1, and if so, we don't care.  We check this at the end of the program.
# Set y=TRUE to keep the outcome; we need this to make sure there is no extra missing data
# in testmodels()
suppressWarnings(
                 samplingmod <- glm(as.formula(samplingformula), x=TRUE, y=TRUE, family=glmlink,
                                    na.action=na.fail, data=data, control=glmcontrol,...)
                 )

# Get the estimated pihats to use for weighting.
data$p.i.h.a.t. <- pihat  <- predict(samplingmod,type='response')

# If any pihats are exactly zero, stop right here -- cannot make estimates!
if ( any(pihat==0) )
  stop(paste(
             "According to the sampling model,",
             "some subjects have zero probability of being observed. \n",
             "Check the sampling model and try again.\n")
       )

# Check if the sampling model converged, and if it didn't, check if it's because some
# pihats are close to 1 or 0.  We expect that if some pihats are close to 1, the routine
# will never claim that it has converged, but it's fine if some pihats are close to 1.
# But if all pihats are not close to either 0 or 1 (close enough as determined by
# criticallogit), then there's a meaningful lack of convergence that we need to warn the
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
# mislead yourself into thinking that the robust variance is the correct variance.
coxmod <- coxph(as.formula(coxformula), weights=1/pihat, method="breslow",
                na.action=na.omit, x=TRUE, subset=TRUE, control=coxphcontrol, data=data,...)


# Get the influence function for the betahats, set to zero for everyone with missing
# covariates. 
D2 <- matrix(0,nrow=n,ncol=length(coxmod$coeff))
D2[observed==1,]  <- residuals(coxmod,'dfbeta',weighted=TRUE)

# Compute the pihat subtract-off by projecting the influences onto the scores.
# Pseudoinverse is important if some combo of J*Delta is always/never observed.
gamma2 <- scores %*% ginv(t(scores)%*%scores) %*% (t(scores) %*% D2)

# This is the correct variance of betahat that accounts for the estimated weights.
D2pihat <- D2-gamma2
varbeta <- t(D2pihat) %*% D2pihat

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

# Invalidate the loglik and score and concordance
coxmod$loglik <- NA
coxmod$score <- NA
coxmod$concordance <- NA

# Plug in the correct Wald test
# This is the usual Wald test, variance estimated at the alternative
# Remove any aliased columns, which will have NA for that coefficient
not.aliased <- which(!is.na(coxmod$coeff))
coeff <- coxmod$coeff[not.aliased]
coxmod$wald.test <- t(coeff) %*% ginv(varbeta[not.aliased,not.aliased]) %*% coeff

##
# Return just cox model or both cox and sampling models
##
if ( outputsamplingmod )
  return(list(coxmod,samplingmod))
else
  return(coxmod)

}
