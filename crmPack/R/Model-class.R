#####################################################################################
## Author: Daniel Sabanes Bove[sabanesd *a*t* roche *.* com],
##         Wai Yin Yeung [ w *.* yeung *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Model-class.R] by DSB Mon 11/05/2015 17:44>
##
## Description:
## Encapsulate the model input in a formal class.
##
## History:
## 31/01/2014   file creation
## 08/07/2015   Adding Pseudo model classes
###################################################################################

##' @include helpers.R
{}


##' Class for All models
##' This is a class where all models inherit.
##' 
##' @slot datanames The names of all data slots that are used in all models. 
##' In particularly, those are also used in the \code{datamodel} and/or
##' \code{priormodel} definition for \code{\linkS4class{GeneralModel}}.
##' 
##' @seealso \code{\linkS4class{GeneralModel}}, \code{\linkS4class{ModelPseudo}}
##' 
##' @export
##' @keywords classes

.AllModels<-setClass(Class="AllModels",
                     representation(datanames="character"))
validObject(.AllModels())

##' No Intitialization function for this

## ============================================================

##' General class for model input
##'
##' This is the general model class, from which all other specific models for 
##' involving BUGS (the software) for computing result. 
##' It inherits all slots from \code{\linkS4class{AllModels}}
##'
##' The \code{datamodel} must obey the convention that the data input is
##' called exactly as in the corresponding data class.
##' All prior distributions for parameters should be contained in the
##' model function \code{priormodel}. The background is that this can
##' be used to simulate from the prior distribution, before obtaining any
##' data.
##'
##' @slot datamodel a function representing the BUGS data model specification
##' (see the details above)
##' @slot priormodel a function representing the BUGS prior specification
##' (see the details above)
##' @slot modelspecs a function computing the list of the data model and prior
##' model specifications that are required for fully specifying them (e.g. prior
##' parameters, reference dose, etc.), based on the data
##' slots that are then required as arguments of this function. This will then
##' be passed to BUGS for the computations.
##' @slot init a function computing the list of starting values for parameters
##' required to be initialized in the MCMC sampler, based on the
##' data slots that are then required as arguments of this
##' function
##' @slot sample names of all parameters from which you would like to save the
##' MCMC samples.
##'
##' @seealso \code{\linkS4class{Model}}
##'
##' @export
##' @keywords classes
.GeneralModel <-
    setClass(Class="GeneralModel",
             representation(datamodel="function",
                            priormodel="function",
                            modelspecs="function",
                            init="function",
                            sample="character"),
             contains="AllModels",
             validity=
               function(object){
                 o <- Validate()
                 
                 o$check(all(names(formals(object@init)) %in%
                               object@datanames),
                         "arguments of the init function must be data names")

                 o$result()
               })
validObject(.GeneralModel())

## no init function for this one



## ============================================================

##' Class for the model input
##'
##' This is the model class for single agent dose escalation,
##' from which all other specific models inherit. It inherits all slots
##' from \code{\linkS4class{GeneralModel}}.
##'
##' The \code{datamodel} must obey the convention that the data input is
##' called exactly as in the \code{\linkS4class{Data}} class.
##' All prior distributions for parameters should be contained in the
##' model function \code{priormodel}. The background is that this can
##' be used to simulate from the prior distribution, before obtaining any
##' data.
##'
##' The \code{dose} function has as first argument \code{prob}, a scalar
##' toxicity probability which is targeted. Additional arguments are model
##' parameters. Then it computes, using model parameter(s) (samples), the
##' resulting dose. Note that the model parameters are called exactly as in the
##' \code{model} and must be included in the \code{sample} vector. The vectors
##' of all samples for these parameters will then be supplied to the function.
##' So your function must be able to process vectors of the model parameters,
##' i.e. it must vectorize over them.
##'
##' The \code{prob} function has as first argument \code{dose}, which is a
##' scalar dose. Additional arguments are model parameters. Then it computes,
##' using model parameter(s) (samples), the resulting probability of toxicity at
##' that dose. Again here, the function must vectorize over the model
##' parameters.
##'
##' If you work with multivariate parameters, then please assume that your
##' the two functions receive either one parameter value as a row vector,
##' or a samples matrix where the rows correspond to the sampling index, i.e.
##' the layout is then nSamples x dimParameter.
##'
##' Note that \code{dose} and \code{prob} are the inverse functions of each
##' other.
##'
##' @slot dose a function computing the dose reaching a specific target
##' probability, based on the model parameters and additional prior settings
##' (see the details above)
##' @slot prob a function computing the probability of toxicity for a specific
##' dose, based on the model parameters and additional prior settings (see the
##' details above)
##'
##' @seealso \code{\linkS4class{LogisticNormal}},
##' \code{\linkS4class{LogisticLogNormal}},
##' \code{\linkS4class{LogisticLogNormalSub}},
##' \code{\linkS4class{LogisticKadane}},
##' \code{\linkS4class{DualEndpoint}}
##'
##' @export
##' @keywords classes
.Model <-
    setClass(Class="Model",
             representation(dose="function",
                            prob="function"),
             contains="GeneralModel",
             validity=
                 function(object){
                     o <- Validate()

                     o$check(all(names(formals(object@dose)) %in%
                                 c("prob", object@sample)),
                             "objects of dose function incorrect")
                     o$check(all(names(formals(object@prob)) %in%
                                 c("dose", object@sample)),
                             "objects of prob function incorrect")

                     o$result()
         })
validObject(.Model())

## no init function for this one

## ============================================================


##' Standard logistic model with bivariate (log) normal prior
##'
##' This is the usual logistic regression model with a bivariate normal prior on
##' the intercept and log slope.
##'
##' The covariate is the natural logarithm of the dose \eqn{x} divided by
##' the reference dose \eqn{x^{*}}:
##'
##' \deqn{logit[p(x)] = \alpha + \beta \cdot \log(x/x^{*})}
##' where \eqn{p(x)} is the probability of observing a DLT for a given dose
##' \eqn{x}.
##'
##' The prior is
##' \deqn{(\alpha, \log(\beta)) \sim Normal(\mu, \Sigma)}
##'
##' The slots of this class contain the mean vector and the covariance matrix of
##' the bivariate normal distribution, as well as the reference dose.
##'
##' Note that the parametrization inside the class uses alpha0 and alpha1.
##' alpha0 is identical to the intercept \eqn{\alpha} above and is the log-odds
##' for a DLT at the reference dose x*. Therefore, the prior mean for alpha0
##' is the expected log-odds at the reference dose x* before observing any data.
##' Note that the expected odds is not just the exp of the prior mean of alpha0,
##' because the non-linearity of the exp transformation. The log-normal
##' distribution on Wikipedia gives the formula for computing the prior mean of
##' exp(alpha0). alpha0 is the log(alpha) in the Neuenschwander et al. (2008)
##' paper. alpha1 is identical to \eqn{\log(\beta)} above and equals the beta
##' in the Neuenschwander et al paper. exp(alpha1) gives the odds-ratio for DLT
##' between two doses that differ by the factor exp(1) ~ 2.7. alpha1 has a
##' log-normal distribution in the LogisticLogNormal model in order to ensure
##' positivity of alpha1 and thus exp(alpha1) > 1.
##'
##' @slot mean the prior mean vector \eqn{\mu}
##' @slot cov the prior covariance matrix \eqn{\Sigma}
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-LogisticLogNormal.R
##' @export
##' @keywords classes
.LogisticLogNormal <-
    setClass(Class="LogisticLogNormal",
             representation(mean="numeric",
                            cov="matrix",
                            refDose="numeric"),
             prototype(mean=c(0, 1),
                       cov=diag(2),
                       refDose=1),
             contains="Model",
             validity=
                 function(object){
                     o <- Validate()

                     o$check(length(object@mean) == 2,
                             "mean must have length 2")
                     o$check(identical(dim(object@cov), c(2L, 2L)) &&
                                 ! is.null(chol(object@cov)),
                             "cov must be positive-definite 2x2 covariance matrix")
                     o$check(is.scalar(object@refDose) &&
                                 (object@refDose > 0),
                             "refDose must be positive scalar")

                     o$result()
                 })
validObject(.LogisticLogNormal())


##' Initialization function for the "LogisticLogNormal" class
##'
##' @param mean the prior mean vector
##' @param cov the prior covariance matrix
##' @param refDose the reference dose
##' @return the \code{\linkS4class{LogisticLogNormal}} object
##'
##' @export
##' @keywords methods
LogisticLogNormal <- function(mean,
                              cov,
                              refDose)
{
    .LogisticLogNormal(mean=mean,
                       cov=cov,
                       refDose=refDose,
                       datamodel=
                       function(){
                           ## the logistic likelihood
                           for (i in 1:nObs)
                           {
                               y[i] ~ dbern(p[i])
                               logit(p[i]) <- alpha0 + alpha1 * StandLogDose[i]
                               StandLogDose[i] <- log(x[i] / refDose)
                           }
                       },
                       priormodel=
                       function(){
                           ## the multivariate normal prior on the (transformed)
                           ## coefficients
                           priorPrec[1:2,1:2] <- inverse(priorCov[,])
                           theta[1:2] ~ dmnorm(priorMean[1:2], priorPrec[1:2,1:2])
                           ## extract actual coefficients
                           alpha0 <- theta[1]
                           alpha1 <- exp(theta[2])

                           ## dummy to use refDose here.
                           ## It is contained in the modelspecs list below,
                           ## so it must occur here
                           bla <- refDose + 1
                       },
                       datanames=c("nObs", "y", "x"),
                       modelspecs=
                       function(){
                           list(refDose=refDose,
                                priorCov=cov,
                                priorMean=mean)
                       },
                       dose=
                       function(prob, alpha0, alpha1){
                           StandLogDose <- (logit(prob) - alpha0) / alpha1
                           return(exp(StandLogDose) * refDose)
                       },
                       prob=
                       function(dose, alpha0, alpha1){
                           StandLogDose <- log(dose / refDose)
                           return(plogis(alpha0 + alpha1 * StandLogDose))
                       },
                       init=
                       ## todo: find better starting values
                       function(){
                           list(theta=c(0, 1))
                       },
                       sample=
                       c("alpha0", "alpha1"))
}
validObject(LogisticLogNormal(mean=c(0, 1),
                              cov=diag(2),
                              refDose=1))


## ============================================================


##' Standard logistic model with bivariate (log) normal prior with substractive
##' dose standardization
##'
##' This is the usual logistic regression model with a bivariate normal prior on
##' the intercept and log slope.
##'
##' The covariate is the dose \eqn{x} minus the reference dose \eqn{x^{*}}:
##'
##' \deqn{logit[p(x)] = \alpha + \beta \cdot (x - x^{*})}
##' where \eqn{p(x)} is the probability of observing a DLT for a given dose
##' \eqn{x}.
##'
##' The prior is
##' \deqn{(\alpha, \log(\beta)) \sim Normal(\mu, \Sigma)}
##'
##' The slots of this class contain the mean vector and the covariance matrix of
##' the bivariate normal distribution, as well as the reference dose.
##'
##' @slot mean the prior mean vector \eqn{\mu}
##' @slot cov the prior covariance matrix \eqn{\Sigma}
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-LogisticLogNormalSub.R
##' @export
##' @keywords classes
.LogisticLogNormalSub <-
    setClass(Class="LogisticLogNormalSub",
             contains="Model",
             representation(mean="numeric",
                            cov="matrix",
                            refDose="numeric"),
             prototype(mean=c(0, 1),
                       cov=diag(2),
                       refDose=1),
             validity=
                 function(object){
                     o <- Validate()

                     o$check(length(object@mean) == 2,
                             "mean must have length 2")
                     o$check(identical(dim(object@cov), c(2L, 2L)) &&
                                 ! is.null(chol(object@cov)),
                             "cov must be positive-definite 2x2 covariance matrix")
                     o$check(is.scalar(object@refDose) &&
                                 (object@refDose > 0),
                             "refDose must be positive scalar")

                     o$result()
                 })
validObject(.LogisticLogNormalSub())


##' Initialization function for the "LogisticLogNormalSub" class
##'
##' @param mean the prior mean vector
##' @param cov the prior covariance matrix
##' @param refDose the reference dose
##' @return the \code{\linkS4class{LogisticLogNormalSub}} object
##'
##' @export
##' @keywords methods
LogisticLogNormalSub <- function (mean,
                                  cov,
                                  refDose){

    .LogisticLogNormalSub(mean=mean,
                          cov=cov,
                          refDose=refDose,
                          datamodel=
                              function(){
                                  ## the logistic likelihood
                                  for (i in 1:nObs)
                                  {
                                      y[i] ~ dbern(p[i])
                                      logit(p[i]) <- alpha0 + alpha1 * StandDose[i]
                                      StandDose[i] <- x[i] - refDose
                                  }
                              },
                          priormodel=
                              function(){
                                  ## the multivariate normal prior on the (transformed)
                                  ## coefficients
                                  priorPrec[1:2,1:2] <- inverse(priorCov[,])
                                  theta[1:2] ~ dmnorm(priorMean[1:2], priorPrec[1:2,1:2])
                                  ## extract actual coefficients
                                  alpha0 <- theta[1]
                                  alpha1 <- exp(theta[2])

                                  ## dummy to use refDose here.
                                  ## It is contained in the modelspecs list below,
                                  ## so it must occur here
                                  bla <- refDose + 1
                              },
                          datanames=c("nObs", "y", "x"),
                          modelspecs=
                              function(){
                                  list(refDose=refDose,
                                       priorCov=cov,
                                       priorMean=mean)
                              },
                          dose=
                              function(prob, alpha0, alpha1){
                                  StandDose <- (logit(prob) - alpha0) / alpha1
                                  return(StandDose + refDose)
                              },
                          prob=
                              function(dose, alpha0, alpha1){
                                  StandDose <- dose - refDose
                                  return(plogis(alpha0 + alpha1 * StandDose))
                              },
                          init=
                              ## todo: find better starting values
                              function(){
                                  list(theta=c(0, -20))
                              },
                          sample=
                              c("alpha0", "alpha1"))
}
validObject(LogisticLogNormalSub(mean=c(0, 1),
                                 cov=diag(2),
                                 refDose=1))

## ============================================================


##' Standard logistic model with bivariate normal prior
##'
##' This is the usual logistic regression model with a bivariate normal prior on
##' the intercept and slope.
##'
##' The covariate is the natural logarithm of the dose \eqn{x} divided by
##' the reference dose \eqn{x^{*}}:
##'
##' \deqn{logit[p(x)] = \alpha + \beta \cdot \log(x/x^{*})}
##' where \eqn{p(x)} is the probability of observing a DLT for a given dose
##' \eqn{x}.
##'
##' The prior is
##' \deqn{(\alpha, \beta) \sim Normal(\mu, \Sigma)}
##'
##' The slots of this class contain the mean vector, the covariance and
##' precision matrices of the bivariate normal distribution, as well as the
##' reference dose.
##'
##' @slot mean the prior mean vector \eqn{\mu}
##' @slot cov the prior covariance matrix \eqn{\Sigma}
##' @slot prec the prior precision matrix \eqn{\Sigma^{-1}}
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-LogisticNormal.R
##' @export
##' @keywords classes
.LogisticNormal <-
    setClass(Class="LogisticNormal",
             contains="Model",
             representation(mean="numeric",
                            cov="matrix",
                            prec="matrix",
                            refDose="numeric"),
             prototype(mean=c(0, 1),
                       cov=diag(2),
                       prec=diag(2),
                       refDose=1),
             validity=
                 function(object){
                     o <- Validate()

                     o$check(length(object@mean) == 2,
                             "mean must have length 2")
                     o$check(identical(dim(object@cov), c(2L, 2L)) &&
                                 ! is.null(chol(object@cov)),
                             "cov must be positive-definite 2x2 covariance matrix")
                     o$check(all.equal(solve(object@cov), object@prec,
                                       check.attributes=FALSE),
                             "prec must be inverse of cov")
                     o$check(is.scalar(object@refDose) &&
                                 (object@refDose > 0),
                             "refDose must be positive scalar")

                     o$result()
                 })


##' Initialization function for the "LogisticNormal" class
##'
##' @param mean the prior mean vector
##' @param cov the prior covariance matrix
##' @param refDose the reference dose
##' @return the \code{\linkS4class{LogisticNormal}} object
##'
##' @export
##' @keywords methods
LogisticNormal <- function (mean,
                            cov,
                            refDose)
{
    .LogisticNormal(mean=mean,
                    cov=cov,
                    prec=solve(cov),
                    refDose=refDose,
                    datamodel=
                        function(){
                            ## the logistic likelihood
                            for (i in 1:nObs)
                            {
                                y[i] ~ dbern(p[i])
                                logit(p[i]) <- alpha0 + alpha1 * StandLogDose[i]
                                StandLogDose[i] <- log(x[i] / refDose)
                            }
                        },
                    priormodel=
                        function(){
                            ## the multivariate normal prior on the coefficients
                            theta[1:2] ~ dmnorm(priorMean[1:2], priorPrec[1:2,1:2])
                            ## extract actual coefficients
                            alpha0 <- theta[1]
                            alpha1 <- theta[2]

                            ## dummy to use refDose here.
                            ## It is contained in the modelspecs list below,
                            ## so it must occur here
                            bla <- refDose + 1
                        },
                    datanames=c("nObs", "y", "x"),
                    modelspecs=
                        function(){
                            list(refDose=refDose,
                                 priorPrec=prec,
                                 priorMean=mean)
                        },
                    dose=
                        function(prob, alpha0, alpha1){
                            StandLogDose <- (logit(prob) - alpha0) / alpha1
                            return(exp(StandLogDose) * refDose)
                        },
                    prob=
                        function(dose, alpha0, alpha1){
                            StandLogDose <- log(dose / refDose)
                                 return(plogis(alpha0 + alpha1 * StandLogDose))
                        },
                    init=
                        ## todo: find better starting values
                        function(){
                            list(theta=c(0, 1))
                        },
                    sample=
                        c("alpha0", "alpha1"))
}
validObject(LogisticNormal(mean=c(0, 1),
                           cov=diag(2),
                           refDose=1))


## ============================================================

##' Reparametrized logistic model
##'
##' This is the logistic model in the parametrization of Kadane et al. (1980).
##'
##' Let \eqn{\rho_{0} = p(x_{min})} be the probability of a DLT and the minimum
##' dose \eqn{x_{min}}, and let \eqn{\gamma} be the dose with target toxicity
##' probability \eqn{\theta}, i.e. \eqn{p(\gamma) = \theta}. Then it can easily
##' be shown that the logistic regression model has intercept
##' \deqn{\frac{\gamma logit(\rho_{0}) - x_{min} logit(\theta)}{\gamma -
##' x_{min}}}{(\gamma logit(\rho_{0}) - x_{min} logit(\theta)) / (\gamma -
##' x_{min})}
##' and slope
##' \deqn{\frac{logit(theta) - logit(\rho_{0})}{\gamma - x_{min}}}{(logit(theta)
##' - logit(\rho_{0})) / (\gamma - x_{min})}
##'
##' The prior is a uniform distribution for \eqn{\gamma} between \eqn{x_{min}}
##' and \eqn{x_{max}}, and for \eqn{\rho_{0}} as well a uniform distribution
##' between \eqn{0} and \eqn{\theta}.
##'
##' The slots of this class, required for creating the model, are the target
##' toxicity, as well as the minimum and maximum of the dose range. Note that
##' these can be different from the minimum and maximum of the dose grid in the
##' data later on.
##'
##' @slot theta the target toxicity probability \eqn{\theta}
##' @slot xmin the minimum of the dose range \eqn{x_{min}}
##' @slot xmax the maximum of the dose range \eqn{x_{max}}
##'
##' @example examples/Model-class-LogisticKadane.R
##' @export
##' @keywords classes
.LogisticKadane <-
    setClass(Class="LogisticKadane",
             contains="Model",
             representation(theta="numeric",
                            xmin="numeric",
                            xmax="numeric"),
             prototype(theta=0.3,
                       xmin=0.1,
                       xmax=1),
             validity=
                 function(object){
                     o <- Validate()

                     o$check(is.probability(object@theta,
                                            bounds=FALSE),
                             "theta must be a probability > 0 and < 1")
                     o$check(object@xmin < object@xmax,
                             "xmin must be smaller than xmax")
                     o$check(is.scalar(object@xmin),
                             "xmin must be scalar")
                     o$check(is.scalar(object@xmax),
                             "xmax must be scalar")

                     o$result()
                 })
validObject(.LogisticKadane())

##' Initialization function for the "LogisticKadane" class
##'
##' @param theta the target toxicity probability
##' @param xmin the minimum of the dose range
##' @param xmax the maximum of the dose range
##' @return the \code{\linkS4class{LogisticKadane}}
##'
##' @export
##' @keywords methods
LogisticKadane <- function(theta,
                           xmin,
                           xmax)
{
    .LogisticKadane(theta=theta,
                    xmin=xmin,
                    xmax=xmax,
                    datamodel=
                        function(){
                            ## the logistic likelihood
                            for (i in 1:nObs)
                            {
                                y[i] ~ dbern(p[i])
                                logit(p[i]) <- (1/(gamma - xmin)) *
                                    (gamma*logit(rho0) - xmin*logit(theta)
                                     + (logit(theta) - logit(rho0)) * x[i])
                            }
                        },
                    priormodel=
                        function(){
                            ## priors
                            gamma ~ dunif(xmin, xmax)
                            rho0 ~ dunif(0, theta)
                        },
                    datanames=c("nObs", "y", "x"),
                    modelspecs=
                        function(){
                            list(theta=theta,
                                 xmin=xmin,
                                 xmax=xmax)
                        },
                    dose=
                        function(prob, rho0, gamma){
                            ret <- gamma * (logit(prob) - logit(rho0)) +
                                xmin * (logit(theta) - logit(prob))
                            ret <- ret / (logit(theta) - logit(rho0))
                            return(ret)
                        },
                    prob=
                        function(dose, rho0, gamma){
                            ret <- (gamma*logit(rho0) - xmin*logit(theta)
                                    + (logit(theta) - logit(rho0)) * dose)
                            ret <- plogis(ret / (gamma - xmin))
                            return(ret)
                        },
                    init=
                        function(){
                            list(rho0 = theta / 10,
                                 gamma = (xmax - xmin) / 2)},
                    sample=
                        c("rho0", "gamma"))
}


## ============================================================


##' Dual endpoint model
##'
##' todo: describe the model
##'
##' @slot mu For the probit toxicity model, \code{mu} contains the prior mean
##' vector
##' @slot Sigma For the probit toxicity model, contains the prior covariance
##' matrix
##' @slot sigma2betaW For the biomarker model, contains the prior variance
##' factor of the random walk prior. If it is not a single number, it can also
##' contain a vector with elements \code{a} and {b} for the inverse-gamma prior
##' on \code{sigma2betaW}.
##' @slot sigma2W Either a fixed value for the biomarker variance, or a vector
##' with elements \code{a} and \code{b} for the inverse-gamma prior parameters.
##' @slot rho Either a fixed value for the correlation (between -1 and 1), or a
##' vector with elements \code{a} and \code{b} for the Beta prior on the
##' transformation kappa = (rho + 1) / 2, which is in (0, 1). For example,
##' \code{a=1,b=1} leads to a uniform prior on rho.
##' @slot useRW1 for specifying the random walk prior on the biomarker level: if
##' \code{TRUE}, RW1 is used, otherwise RW2.
##' @slot useFixed a list with logical value for each of the three parameters
##' \code{sigma2betaW}, \code{sigma2W} and \code{rho} indicating whether
##' a fixed value is used or not.
##'
##' @export
##' @keywords classes internal
setClass(Class="DualEndpointOld",
         contains="Model",
         representation=
         representation(mu="numeric",
                        Sigma="matrix",
                        sigma2betaW="numeric",
                        sigma2W="numeric",
                        rho="numeric",
                        useRW1="logical",
                        useFixed="list"),
         validity=
         function(object){

             ## check the prior parameters with variable content
             for(parName in c("sigma2betaW", "sigma2W", "rho"))
             {
                 ## if we use a fixed value for this parameter
                 if(object@useFixed[[parName]])
                 {
                     ## check range of value
                     if(parName == "rho")
                     {
                         stopifnot((object@rho > -1) && (object@rho < 1))
                     } else {
                         stopifnot(slot(object, parName) > 0)
                     }
                 } else {
                     ## use a IG(a, b) or Beta(a, b)  prior
                     stopifnot(all(slot(object, parName) > 0),
                               identical(names(slot(object, parName)),
                                         c("a", "b")))
                 }
             }

             ## check the other prior parameters
             stopifnot(identical(length(object@mu), 2L),
                       identical(dim(object@Sigma), c(2L, 2L)),
                       is.scalar(object@useRW1))
         })


##' Initialization method for the "DualEndpointOld" class
##'
##' @param .Object the \code{\linkS4class{DualEndpointOld}} we want to
##' initialize
##' @param mu see \code{\linkS4class{DualEndpointOld}}
##' @param Sigma see \code{\linkS4class{DualEndpointOld}}
##' @param sigma2betaW see \code{\linkS4class{DualEndpointOld}}
##' @param sigma2W see \code{\linkS4class{DualEndpointOld}}
##' @param rho see \code{\linkS4class{DualEndpointOld}}
##' @param smooth either \dQuote{RW1} (default) or \dQuote{RW2}, for
##' specifying the random walk prior on the biomarker level.
##' @param \dots not used
##'
##' @export
##' @keywords methods
setMethod("initialize",
          signature(.Object = "DualEndpointOld"),
          function(.Object,
                   mu,
                   Sigma,
                   sigma2betaW,
                   sigma2W,
                   rho,
                   smooth=c("RW1", "RW2"),
                   ...){

              ## Find out RW choice
              smooth <- match.arg(smooth)
              .Object@useRW1 <- smooth == "RW1"

              ## Find out which parameters are fixed
              useFixed <- list()
              for(parName in c("sigma2betaW", "sigma2W", "rho"))
              {
                  useFixed[[parName]] <-
                      identical(length(get(parName)), 1L)
              }
              .Object@useFixed <- useFixed

              ## build together the prior model and the parameters
              ## to be saved during sampling
              ## ----------

              ## start with this:

              modelspecs <-
                  list(mu=mu,
                       PrecBetaZ=solve(Sigma),
                       low=c(-10000, 0),
                       high=c(0, 10000))

              priormodel <-
                  function(){
                      ## priors
                      betaW[1:nGrid] ~ car.normal(adj[], weights[],
                                                  num[], precBetaW)
                      ## note that "adj", "weights" and "num" come
                      ## from the data object! Because it depends on
                      ## the dose grid
                      betaWintercept ~ dflat()

                      ## the bivariate normal prior for the
                      ## probit coefficients
                      betaZ[1:2] ~ dmnorm(mu[], PrecBetaZ[,])

                      ## conditional precision for biomarker
                      condPrecW <- precW / (1 - pow(rho, 2))
                  }

              ## we will fill in more, depending on which parameters
              ## are fixed, in these two variables:
              sample <- c("betaZ", "betaWintercept", "betaW")
              initlist <- list()

              ## first the biomarker regression variance
              if(! useFixed[["sigma2W"]])
              {
                  priormodel <-
                      joinModels(priormodel,
                                 function(){
                                     ## gamma prior for biomarker precision
                                     precW ~ dgamma(precWa, precWb)
                                 })

                  sample <- c(sample,
                              "precW")

                  modelspecs <- c(modelspecs,
                                  list(precWa=sigma2W["a"],
                                       precWb=sigma2W["b"]))

                  initlist$precW <- 1

              } else {
                  modelspecs <- c(modelspecs,
                                  list(precW=1/sigma2W))
              }

              ## second the variance for the RW prior
              if(! useFixed[["sigma2betaW"]])
              {
                  priormodel <-
                      joinModels(priormodel,
                                 function(){
                                     ## gamma prior for RW precision
                                     precBetaW ~ dgamma(precBetaWa, precBetaWb)
                                 })

                  sample <- c(sample,
                              "precBetaW")

                  initlist$precBetaW <- 1

                  modelspecs <- c(modelspecs,
                                  list(precBetaWa=sigma2betaW["a"],
                                       precBetaWb=sigma2betaW["b"]))
              } else {
                  modelspecs <- c(modelspecs,
                                  list(precBetaW=1/sigma2betaW))
              }

              ## third the correlation
              if(! useFixed[["rho"]])
              {
                  priormodel <-
                      joinModels(priormodel,
                                 function(){
                                     ## transformed Beta prior for rho
                                     kappa ~ dbeta(rhoa, rhob)
                                     rho <- kappa * 2 - 1
                                 })

                  sample <- c(sample,
                              "rho")

                  initlist$kappa <- 1/2

                  modelspecs <- c(modelspecs,
                                  list(rhoa=rho["a"],
                                       rhob=rho["b"]))
              } else {
                  modelspecs <- c(modelspecs,
                                  list(rho=rho))
              }

              ## now build the *function* modelspecs, which computes
              ## from data slots the required RW matrices etc.
              modelspecsFun <- function(x, xLevel, nGrid)
              {
                  ## design matrices for tox and biomarker models:

                  ## tox
                  designZ <- cbind(1, x)
                  ## todo: note that this is the easiest case here.
                  ## we could in principle employ any probit regression model
                  ## for the toxicity! So later on, we can extend this
                  ## to make it more flexible.

                  ## biomarker
                  designW <-
                      model.matrix(~ - 1 +
                                   I(factor(xLevel,
                                            levels=seq_len(nGrid))))
                  dimnames(designW) <- list(NULL, NULL)

                  ## difference matrix of order 1:
                  D1mat <- cbind(0,
                                 diag(nGrid - 1)) -
                                     cbind(diag(nGrid - 1), 0)

                  ## we will compute the RW prior data in the following
                  weights <- adj <- num <- numeric()
                  ## Note that these vectors are sized during the computations
                  ## now.

                  if(.Object@useRW1)
                  {
                      ## RW1
                      ## ----------

                      ## the rank-deficient prior precision for the RW1 prior:
                      RWmat <- crossprod(D1mat)
                      ## Note that this needs to be divided by sigma2betaW to
                      ## obtain final prior precision

                      ## the rank of this matrix
                      RWmatRank <- nGrid - 1L

                      ## compute vectors
                      for(t in 1:1)
                      {
                          weights[t] <- 1;
                          adj[t] <- t+1;
                          num[t] <- 1
                      }
                      for(t in 2:(nGrid-1))
                      {
                          weights[2+(t-2)*2] <- 1;
                          adj[2+(t-2)*2] <- t-1
                          weights[3+(t-2)*2] <- 1;
                          adj[3+(t-2)*2] <- t+1;
                          num[t] <- 2
                      }
                      for(t in nGrid:nGrid)
                      {
                          weights[(nGrid-2)*2 + 2] <- 1;
                          adj[(nGrid-2)*2 + 2] <- t-1;
                          num[t] <- 1
                      }

                  } else {

                      ## RW2
                      ## ----------

                      ## for second-order differences:
                      D2mat <- D1mat[-1, -1] %*% D1mat

                      ## same for RW2
                      RWmat <- crossprod(D2mat)
                      RWmatRank <- nGrid - 2L

                      ## compute vectors
                      for(t in 1:1) {
                          weights[t] <- 2; adj[t] <- t+1
                          weights[t+1] <- -1; adj[t+1] <- t+2; num[t] <- 2
                      }
                      for(t in 2:2) {
                          weights[t+1] <- 2; adj[t+1] <- t-1
                          weights[t+2] <- 4; adj[t+2] <- t+1
                          weights[t+3] <- -1; adj[t+3] <- t+2; num[t] <- 3
                      }
                      for(t in 3:(nGrid-2)) {
                          weights[6+(t-3)*4] <- -1; adj[6+(t-3)*4] <- t-2
                          weights[7+(t-3)*4] <- 4; adj[7+(t-3)*4] <- t-1
                          weights[8+(t-3)*4] <- 4; adj[8+(t-3)*4] <- t+1
                          weights[9+(t-3)*4] <- -1; adj[9+(t-3)*4] <- t+2;
                          num[t] <- 4
                      }
                      for(t in (nGrid-1):(nGrid-1)) {
                          weights[(nGrid-4)*4 + 6] <- 2;
                          adj[(nGrid-4)*4 + 6] <- t+1
                          weights[(nGrid-4)*4 + 7] <- 4;
                          adj[(nGrid-4)*4 + 7] <- t-1
                          weights[(nGrid-4)*4 + 8] <- -1;
                          adj[(nGrid-4)*4 + 8] <- t-2;
                          num[t] <- 3
                      }
                      for(t in nGrid:nGrid) {
                          weights[(nGrid-4)*4 + 9] <- 2;
                          adj[(nGrid-4)*4 + 9] <- t-1
                          weights[(nGrid-4)*4 + 10] <- -1;
                          adj[(nGrid-4)*4 + 10] <- t-2;
                          num[t] <- 2
                      }
                  }

                  ## finally return the list
                  return(c(modelspecs,
                           list(## designZ=designZ,
                                ## designW=designW,
                                ## RWmat=RWmat,
                                ## RWmatRank=RWmatRank,
                                weights=weights,
                                adj=adj,
                                num=num)))
              }


              ## go to the general initialize method now
              callNextMethod(.Object,
                             mu=mu,
                             Sigma=Sigma,
                             sigma2betaW=sigma2betaW,
                             sigma2W=sigma2W,
                             rho=rho,
                             datamodel=
                             function(){
                                 ## the likelihood
                                 for (i in 1:nObs)
                                 {
                                     ## the toxicity model
                                     z[i] ~ dnorm(meanZ[i], 1) %_%
                                         I(low[y[i] + 1], high[y[i] + 1])

                                     ## the conditional biomarker model
                                     w[i] ~ dnorm(condMeanW[i], condPrecW)

                                     ## the moments
                                     meanZ[i] <- betaZ[1] + betaZ[2] * x[i]
                                     condMeanW[i] <- betaWintercept + betaW[xLevel[i]] +
                                         rho / sqrt(precW) * (z[i] - meanZ[i])
                                     ## Note that betaW has a sum to zero constraint here.
                                     ## Therefore we have to add an intercept with a flat prior
                                     ## on top.
                                 }
                             },
                             priormodel=priormodel,
                             datanames=
                             c("nObs", "w", "x", "xLevel", "y", "nGrid"),
                             modelspecs=modelspecsFun,
                             dose=
                             function(prob, betaZ){
                                 ret <- (qnorm(prob) - betaZ[, 1]) / betaZ[, 2]
                                 return(ret)
                             },
                             prob=
                             function(dose, betaZ){
                                 ret <- pnorm(betaZ[, 1] + betaZ[, 2] * dose)
                                 return(ret)
                             },
                             init=
                             function(y, w, nGrid){
                                 c(initlist,
                                   list(z=
                                        ifelse(y==0, -1, 1),
                                        betaZ=c(0,1),
                                        betaWintercept=mean(w),
                                        betaW=
                                        rep(0, nGrid)))},
                             sample=sample,
                             ...)
          })


## ============================================================


##' General class for the dual endpoint model
##'
##' The idea of the dual-endpoint models is to model not only the dose-toxicity
##' relationship, but also to model at the same time the relationship of a PD
##' biomarker with the dose. The subclasses of this class detail how the
##' dose-biomarker relationship is parametrized and are those to be used. This
##' class here shall contain all the common features to reduce duplicate code.
##' (However, this class must not be virtual, because we need to create objects
##' of it during the construction of subclass objects.)
##'
##' Currently a probit regression model
##' \deqn{\Phi^{-1}[p(x)] = \beta_{Z1} + \beta_{Z2} \cdot x}
##' is used, where \eqn{p(x)} is the probability of observing a DLT for a given
##' dose \eqn{x}, and \eqn{\Phi} is the standard normal cdf. This could later be
##' generalized to have a reference dose or a log transformation for the dose.
##' The prior is
##' \deqn{\left( \beta_{Z1} , log(\beta_{Z2}) \right) \sim Normal(\mu, \Sigma)}.
##'
##' For the biomarker response w at a dose x, we assume
##' \deqn{w(x) \sim Normal(f(x), \sigma^{2}_{W})}
##' and \eqn{f(x)} is a function of the dose x, which is further specified in
##' the subclasses. The biomarker variance \eqn{\sigma^{2}_{W}} can be fixed or
##' assigned an inverse gamma prior distribution; see the details below under
##' slot \code{sigma2W}.
##'
##' Finally, the two endpoints y (the binary DLT variable) and w (the biomarker)
##' can be correlated, by assuming a correlation \eqn{\rho} between the
##' underlying continuous latent toxicity variable z and the biomarker w.
##' Again, this correlation can be fixed or assigned a prior distribution from
##' the scaled beta family; see the details below under slot \code{rho}.
##'
##' Please see the Hive page for more details on the model and the example
##' vignette by typing \code{crmPackExample()} for a full example.
##'
##' @slot mu For the probit toxicity model, \code{mu} contains the prior mean
##' vector
##' @slot Sigma For the probit toxicity model, contains the prior covariance
##' matrix
##' @slot sigma2W Either a fixed value for the biomarker variance, or a vector
##' with elements \code{a} and \code{b} for the inverse-gamma prior parameters.
##' @slot rho Either a fixed value for the correlation (between -1 and 1), or a
##' vector with elements \code{a} and \code{b} for the Beta prior on the
##' transformation kappa = (rho + 1) / 2, which is in (0, 1). For example,
##' \code{a=1,b=1} leads to a uniform prior on rho.
##' @slot useFixed a list with logical value for each of the parameters
##' indicating whether a fixed value is used or not; this slot is needed for
##' internal purposes and not to be touched by the user.
##'
##' @export
##' @seealso Current subclasses: \code{\linkS4class{DualEndpointRW}},
##' \code{\linkS4class{DualEndpointBeta}}
##' @keywords classes
.DualEndpoint <-
    setClass(Class="DualEndpoint",
             representation(mu="numeric",
                            Sigma="matrix",
                            sigma2W="numeric",
                            rho="numeric",
                            useFixed="list"),
             prototype(mu=c(0, 1),
                       Sigma=diag(2),
                       sigma2W=1,
                       rho=0,
                       useFixed=
                           list(sigma2W=TRUE,
                                rho=TRUE)),
             contains="Model",
             validity=
                 function(object){
                     o <- Validate()

                     ## check the prior parameters with variable content
                     for(parName in c("sigma2W", "rho"))
                     {
                         ## if we use a fixed value for this parameter
                         if(object@useFixed[[parName]])
                         {
                             ## check range of value
                             if(parName == "rho")
                             {
                                 o$check((object@rho > -1) && (object@rho < 1),
                                         "rho must be in (-1, 1)")
                             } else {
                                 o$check(slot(object, parName) > 0,
                                         paste(parName, "must be positive"))
                             }
                         } else {
                             ## use a IG(a, b) or Beta(a, b)  prior
                             o$check(identical(names(slot(object, parName)),
                                               c("a", "b")),
                                     paste(parName,
                                           "must have names 'a' and 'b'"))
                             o$check(all(slot(object, parName) > 0),
                                     paste(parName,
                                           "must have positive prior parameters"))
                         }
                     }

                     ## check the other prior parameters
                     o$check(identical(length(object@mu), 2L),
                             "mu must have length 2")
                     o$check(identical(dim(object@Sigma), c(2L, 2L)) &&
                                 ! is.null(chol(object@Sigma)),
                             "Sigma must be positive-definite 2x2 covariance matrix")

                     o$result()
                 })
validObject(.DualEndpoint())

##' Initialization function for the "DualEndpoint" class
##'
##' @param mu see \code{\linkS4class{DualEndpoint}}
##' @param Sigma see \code{\linkS4class{DualEndpoint}}
##' @param sigma2W see \code{\linkS4class{DualEndpoint}}
##' @param rho see \code{\linkS4class{DualEndpoint}}
##' @return the \code{\linkS4class{DualEndpoint}} object
##'
##' @export
##' @keywords methods
DualEndpoint <- function(mu,
                         Sigma,
                         sigma2W,
                         rho)
{
    ## Find out which parameters are fixed
    useFixed <- list()
    for(parName in c("sigma2W", "rho"))
    {
        useFixed[[parName]] <-
            identical(length(get(parName)), 1L)
    }

    ## build together the prior model and the parameters
    ## to be saved during sampling
    ## ----------

    ## start with this:

    modelspecs <-
        list(mu=mu,
             PrecBetaZ=solve(Sigma)## ,
             ## low=c(-10000, 0),
             ## high=c(0, 10000)
             )

    priormodel <-
        function(){
            ## priors for betaW: defined in subclasses

            ## the bivariate normal prior for the
            ## probit coefficients
            log.betaZ[1:2] ~ dmnorm(mu[], PrecBetaZ[,])
            
            betaZ[1] <- log.betaZ[1]
            betaZ[2] <- exp(log.betaZ[2])

            ## conditional precision for biomarker
            condPrecW <- precW / (1 - pow(rho, 2))
        }

    ## we will fill in more, depending on which parameters
    ## are fixed, in these two variables:
    sample <- c("betaZ")
    initlist <- list()

    ## first the biomarker regression variance
    if(! useFixed[["sigma2W"]])
    {
        priormodel <-
            joinModels(priormodel,
                       function(){
                           ## gamma prior for biomarker precision
                           precW ~ dgamma(precWa, precWb)
                       })

        sample <- c(sample,
                    "precW")

        modelspecs <- c(modelspecs,
                        list(precWa=sigma2W["a"],
                             precWb=sigma2W["b"]))

        initlist$precW <- 1

    } else {
        modelspecs <- c(modelspecs,
                        list(precW=1/sigma2W))
    }

    ## second the correlation
    if(! useFixed[["rho"]])
    {
        priormodel <-
            joinModels(priormodel,
                       function(){
                           ## transformed Beta prior for rho
                           kappa ~ dbeta(rhoa, rhob)
                           rho <- kappa * 2 - 1
                       })

        sample <- c(sample,
                    "rho")

        initlist$kappa <- 1/2

        modelspecs <- c(modelspecs,
                        list(rhoa=rho["a"],
                             rhob=rho["b"]))
    } else {
        modelspecs <- c(modelspecs,
                        list(rho=rho))
    }

    ## finally build the object
    .DualEndpoint(mu=mu,
                  Sigma=Sigma,
                  sigma2W=sigma2W,
                  rho=rho,
                  useFixed=useFixed,
                  datamodel=
                  function(){
                      ## the likelihood
                      for (i in 1:nObs)
                      {
                          ## the toxicity model
                          ## z[i] ~ dnorm(meanZ[i], 1) %_%
                          ##     I(low[y[i] + 1], high[y[i] + 1])
                          y[i] ~ dinterval(z[i], 0)
                          z[i] ~ dnorm(meanZ[i], 1)

                          ## the conditional biomarker model
                          w[i] ~ dnorm(condMeanW[i], condPrecW)

                          ## the moments
                          meanZ[i] <- betaZ[1] + betaZ[2] * x[i]
                          condMeanW[i] <- betaW[xLevel[i]] +
                              rho / sqrt(precW) * (z[i] - meanZ[i])
                          ## betaW needs to be defined in subclasses!
                      }
                  },
                  priormodel=priormodel,
                  datanames=
                  c("nObs", "w", "x", "xLevel", "y", "nGrid"),
                  modelspecs=
                  function(){
                      modelspecs
                  },
                  dose=
                  function(prob, betaZ){
                      ret <- (qnorm(prob) - betaZ[, 1]) / betaZ[, 2]
                      return(ret)
                  },
                  prob=
                  function(dose, betaZ){
                      ret <- pnorm(betaZ[, 1] + betaZ[, 2] * dose)
                      return(ret)
                  },
                  init=
                  function(y, w, nGrid){
                      c(initlist,
                        list(z=
                             ifelse(y==0, -1, 1),
                             log.betaZ=c(0,1)))},
                  sample=sample)
}
validObject(DualEndpoint(mu=c(0, 1),
                         Sigma=diag(2),
                         sigma2W=1,
                         rho=0))


## ============================================================


##' Dual endpoint model with RW prior for biomarker
##'
##' This class extends the \code{\linkS4class{DualEndpoint}} class. Here the
##' dose-biomarker relationship \eqn{f(x)} is modelled by a non-parametric
##' random-walk of first (RW1) or second order (RW2) (todo: warning: at the
##' moment only the first order random walk produces useful results).
##'
##' That means, for the RW1 we assume
##' \deqn{\beta_{W,i} - \beta_{W,i-1} \sim Normal(0, \sigma^{2}_{\beta_{W}}),}
##' where \eqn{\beta_{W,i} = f(x_{i})} is the biomarker mean at the i-th dose
##' gridpoint \eqn{x_{i}}.
##' For the RW2, the second-order differences instead of the first-order
##' differences of the biomarker means follow the normal distribution.
##'
##' The variance parameter \eqn{\sigma^{2}_{\beta_{W}}} is important because it
##' steers the smoothness of the function f(x): if it is large, then f(x) will
##' be very wiggly; if it is small, then f(x) will be smooth. This parameter can
##' either be fixed or assigned an inverse gamma prior distribution.
##'
##' Usually this modelling will only make sense if a regular dose grid is used,
##' with equidistant grid points ensuring that the distance \eqn{x_{i} -
##' x_{i-1}} is the same for all grid positions \eqn{i}.
##'
##' @slot sigma2betaW Contains the prior variance factor of the random walk
##' prior for the biomarker model. If it is not a single number, it can also
##' contain a vector with elements \code{a} and {b} for the inverse-gamma prior
##' on \code{sigma2betaW}.
##' @slot useRW1 for specifying the random walk prior on the biomarker level: if
##' \code{TRUE}, RW1 is used, otherwise RW2.
##'
##' @example examples/Model-class-DualEndpointRW.R
##' @export
##' @keywords classes
.DualEndpointRW <-
    setClass(Class="DualEndpointRW",
             representation(sigma2betaW="numeric",
                            useRW1="logical"),
             prototype(sigma2betaW=1,
                       useRW1=TRUE,
                       useFixed=
                       list(sigma2W=TRUE,
                            rho=TRUE,
                            sigma2betaW=TRUE)),
             contains="DualEndpoint",
             validity=
                 function(object){
                     o <- Validate()

                     ## check the additional prior parameters
                     for(parName in c("sigma2betaW"))
                     {
                         ## if we use a fixed value for this parameter
                         if(object@useFixed[[parName]])
                         {
                             o$check(slot(object, parName) > 0,
                                     paste(parName, "must be positive"))
                         } else {
                              o$check(identical(names(slot(object, parName)),
                                               c("a", "b")),
                                     paste(parName,
                                           "must have names 'a' and 'b'"))
                              o$check(all(slot(object, parName) > 0),
                                      paste(parName,
                                            "must have positive prior parameters"))
                         }
                     }

                     o$result()
                 })
validObject(.DualEndpointRW())


##' Initialization function for the "DualEndpointRW" class
##'
##' @param sigma2betaW see \code{\linkS4class{DualEndpointRW}}
##' @param smooth either \dQuote{RW1} (default) or \dQuote{RW2}, for
##' specifying the random walk prior on the biomarker level.
##' @param \dots additional parameters, see \code{\linkS4class{DualEndpoint}}
##' @return the \code{\linkS4class{DualEndpointRW}} object
##'
##' @export
##' @keywords methods
DualEndpointRW <- function(sigma2betaW,
                           smooth=c("RW1", "RW2"),
                           ...)
{
    ## call the initialize function from DualEndpoint
    ## to get started
    start <- DualEndpoint(...)

    ## Find out RW choice
    smooth <- match.arg(smooth)
    useRW1 <- smooth == "RW1"

    if(! useRW1)
    {
        warning("todo: Currently only the RW1 model produces useful results!")
    }

    ## Find out which of the additional parameters are fixed
    for(parName in c("sigma2betaW"))
    {
        start@useFixed[[parName]] <-
            identical(length(get(parName)), 1L)
    }

    ## build together the prior model and the parameters
    ## to be saved during sampling
    ## ----------

    start@priormodel <-
        joinModels(start@priormodel,
                   function(){

                       betaW[1] <- betaWintercept
                       for (j in 2:nGrid) {
                           betaW[j] <- betaWintercept +
                               sum(delta[1:(j-1)])
                       }

                       ## delta will then be defined below
                       ## (depending on whether RW1 or RW2 is used)

                       ## the intercept (= first location betaW value)
                       betaWintercept ~ dnorm(0, 0.000001)
                       ## ~ essentially dflat(),
                       ## which is not available in JAGS.
                   })

    if(useRW1)
    {
        ## add RW1 part
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## the iid first oder differences:
                           for (j in 2:nGrid) {
                               delta[j-1] ~ dnorm(0, precBetaW)
                           }
                       })
    } else {
        ## add RW2 part
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## the first order differences:
                           delta[1] <- deltaStart
                           for (j in 2:(nGrid-1)) {
                               delta[j] <- deltaStart +
                                   sum(delta2[1:(j-1)])
                           }

                           ## the iid second oder differences:
                           for (j in 1:(nGrid-2)) {
                               delta2[j] ~ dnorm(0, precBetaW)
                           }

                           ## the first 1st order difference:
                           deltaStart ~ dnorm(0, 0.000001)
                       })
    }

    ## we will fill in more, depending on which parameters
    ## are fixed, in these two variables:
    start@sample <- c(start@sample,
                      "betaWintercept", "betaW", "delta")

    ## copy some things that we need to amend
    oldModelspecs <- start@modelspecs
    oldInit <- start@init

    ## check the variance for the RW prior
    if(! start@useFixed[["sigma2betaW"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## gamma prior for RW precision
                           precBetaW ~ dgamma(precBetaWa, precBetaWb)
                       })

        start@sample <- c(start@sample,
                          "precBetaW")

        start@init <-
            function(y, w, nGrid){
                c(oldInit(y, w, nGrid),
                  list(precBetaW=1))
            }

        start@modelspecs <-
            function(){
                c(oldModelspecs(),
                  list(precBetaWa=sigma2betaW["a"],
                       precBetaWb=sigma2betaW["b"]))
            }

    } else {
        start@modelspecs <-
            function(){
                c(oldModelspecs(),
                  list(precBetaW=1/sigma2betaW))
            }
    }

    ## call the constructor:
    ## inherit everything from DualEndpoint object "start", and add
    ## specifics
    .DualEndpointRW(start,
                    sigma2betaW=sigma2betaW,
                    useRW1=useRW1)
}
validObject(DualEndpointRW(sigma2betaW=1,
                           smooth="RW1",
                           mu=c(0, 1),
                           Sigma=diag(2),
                           sigma2W=1,
                           rho=0))


## ============================================================


##' Dual endpoint model with beta function for dose-biomarker relationship
##'
##' This class extends the \code{\linkS4class{DualEndpoint}} class. Here the
##' dose-biomarker relationship \eqn{f(x)} is modelled by a parametric, rescaled
##' beta density function:
##'
##' \deqn{f(x) = E_{0} + (E_{max} - E_{0}) * Beta(\delta_{1}, \delta_{2}) *
##'                   (x/x^{*})^{\delta_{1}} * (1 - x/x^{*})^{\delta_{2}}}
##'
##' where \eqn{x^{*}} is the maximum dose (end of the dose range to be
##' considered), \eqn{\delta_{1}} and \eqn{\delta_{2}} are the two beta
##' parameters, and \eqn{E_{0}} and \eqn{E_{max}} are the minimum and maximum
##' levels, respectively. For ease of interpretation, we parametrize with
##' \eqn{\delta_{1}} and the mode of the curve instead, where
##' \deqn{mode = \delta_{1} / (\delta_{1} + \delta_{2}),}
##' and multiplying this with \eqn{x^{*}} gives the mode on the dose grid.
##'
##' All parameters can currently be assigned uniform distributions or be fixed
##' in advance.
##'
##' @slot E0 either a fixed number or the two uniform distribution parameters
##' @slot Emax either a fixed number or the two uniform distribution parameters
##' @slot delta1 either a fixed number or the two uniform distribution parameters
##' @slot mode either a fixed number or the two uniform distribution parameters
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-DualEndpointBeta.R
##' @export
##' @keywords classes
.DualEndpointBeta <-
    setClass("DualEndpointBeta",
             representation(E0="numeric",
                            Emax="numeric",
                            delta1="numeric",
                            mode="numeric",
                            refDose="numeric"),
             prototype(E0=c(0, 100),
                       Emax=c(0, 500),
                       delta1=c(0, 5),
                       mode=c(1, 15),
                       refDose=1000,
                       useFixed=
                       list(sigma2W=TRUE,
                            rho=TRUE,
                            E0=FALSE,
                            Emax=FALSE,
                            delta1=FALSE,
                            mode=FALSE)),
             contains="DualEndpoint",
             validity=
                 function(object){
                     o <- Validate()

                     ## check the prior parameters with variable content
                     for(parName in c("E0", "Emax", "delta1", "mode"))
                     {
                         ## if we use a fixed value for this parameter
                         if(object@useFixed[[parName]])
                         {
                             ## check range of value
                             o$check(slot(object, parName) > 0,
                                         paste(parName, "must be positive"))
                         } else {
                             ## use a Uniform(a, b) prior
                             o$check(all(slot(object, parName) >= 0) &&
                                         (diff(slot(object, parName)) > 0),
                                     paste(parName,
                                           "has not proper prior parameters"))
                         }
                     }

                     ## check the refDose
                     o$check(object@refDose > 0,
                             "refDose must be positive")

                     o$result()
                 })
validObject(.DualEndpointBeta())

##' Initialization function for the "DualEndpointBeta" class
##'
##' @param E0 see \code{\linkS4class{DualEndpointBeta}}
##' @param Emax see \code{\linkS4class{DualEndpointBeta}}
##' @param delta1 see \code{\linkS4class{DualEndpointBeta}}
##' @param mode see \code{\linkS4class{DualEndpointBeta}}
##' @param refDose see \code{\linkS4class{DualEndpointBeta}}
##' @param \dots additional parameters, see \code{\linkS4class{DualEndpoint}}
##' @return the \code{\linkS4class{DualEndpointBeta}} object
##'
##' @export
##' @keywords methods
DualEndpointBeta <- function(E0,
                             Emax,
                             delta1,
                             mode,
                             refDose,
                             ...)
{
    ## call the initialize function from DualEndpoint
    ## to get started
    start <- DualEndpoint(...)

    ## we need the dose grid here in the BUGS model,
    ## therefore add it to datanames
    start@datanames <- c(start@datanames,
                         "doseGrid")

    ## Find out which of the additional parameters are fixed
    for(parName in c("E0", "Emax", "delta1", "mode"))
    {
        start@useFixed[[parName]] <-
            identical(length(get(parName)), 1L)
    }

    ## build together the prior model and the parameters
    ## to be saved during sampling
    ## ----------

    start@priormodel <-
        joinModels(start@priormodel,
                   function(){
                       ## delta2 <- delta1 * (1 - (mode/refDose)) / (mode/refDose)
                       delta2 <- delta1 * (refDose/mode - 1)
                       ## betafun <- (delta1 + delta2)^(delta1 + delta2) *
                       ##     delta1^(- delta1) * delta2^(- delta2)
                       betafun <- (1 + delta2 / delta1)^delta1 *
                           (delta1 / delta2 + 1)^delta2

                       for (j in 1:nGrid)
                       {
                           StandDose[j] <- doseGrid[j] / refDose
                           betaW[j] <- E0 + (Emax - E0) * betafun *
                               StandDose[j]^delta1 * (1 - StandDose[j])^delta2
                       }
                   })

    ## we will fill in more, depending on which parameters
    ## are fixed, in these two variables:
    start@sample <- c(start@sample,
                      "betaW")
    newInits <- list()
    newModelspecs <- list(refDose=refDose)

    ## for E0:
    if(! start@useFixed[["E0"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for E0
                           E0 ~ dunif(E0low, E0high)
                       })

        start@sample <- c(start@sample,
                          "E0")

        newInits$E0 <- mean(E0)
        newModelspecs$E0low <- E0[1]
        newModelspecs$E0high <- E0[2]
    } else {
        newModelspecs$E0 <- E0
    }

    ## for Emax:
    if(! start@useFixed[["Emax"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for Emax
                           Emax ~ dunif(EmaxLow, EmaxHigh)
                       })

        start@sample <- c(start@sample,
                          "Emax")

        newInits$Emax <- mean(Emax)
        newModelspecs$EmaxLow <- Emax[1]
        newModelspecs$EmaxHigh <- Emax[2]
    } else {
        newModelspecs$Emax <- Emax
    }

    ## for delta1 and delta2:
    if(! start@useFixed[["delta1"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for E0
                           delta1 ~ dunif(delta1Low, delta1High)
                       })

        start@sample <- c(start@sample,
                          "delta1")

        newInits$delta1 <- mean(delta1)
        newModelspecs$delta1Low <- delta1[1]
        newModelspecs$delta1High <- delta1[2]
    } else {
        newModelspecs$delta1 <- delta1
    }

    if(! start@useFixed[["mode"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for E0
                           mode ~ dunif(modeLow, modeHigh)
                       })

        start@sample <- c(start@sample,
                          "mode")

        newInits$mode <- mean(mode)
        newModelspecs$modeLow <- mode[1]
        newModelspecs$modeHigh <- mode[2]
    } else {
        newModelspecs$mode <- mode
    }

    ## now define the new modelspecs and init functions:
    oldModelspecs <- start@modelspecs
    start@modelspecs <- function()
    {
        c(oldModelspecs(),
          newModelspecs)
    }

    oldInit <- start@init
    start@init <- function(y, w, nGrid)
    {
        c(oldInit(y, w, nGrid),
          newInits)
    }

    ## finally call the constructor
    .DualEndpointBeta(start,
                      E0=E0,
                      Emax=Emax,
                      delta1=delta1,
                      mode=mode,
                      refDose=refDose)
}
validObject(DualEndpointBeta(E0=10,
                             Emax=50,
                             delta1=c(1, 5),
                             mode=c(3, 10),
                             refDose=10,
                             mu=c(0, 1),
                             Sigma=diag(2),
                             sigma2W=1,
                             rho=0))


## ============================================================


##' Dual endpoint model with emax function for dose-biomarker relationship
##'
##' This class extends the \code{\linkS4class{DualEndpoint}} class. Here the
##' dose-biomarker relationship \eqn{f(x)} is modelled by a parametric EMAX function:
##'
##' \deqn{f(x) = E_{0} + \frac{(E_{max} - E_{0}) * (x/x^{*})}{ED_{50} + (x/x^{*})}}
##'
##' where \eqn{x^{*}} is a reference dose, \eqn{E_{0}} and \eqn{E_{max}} are the
##' minimum and maximum levels for the biomarker and \eqn{ED_{50}} is the dose
##' achieving half of the maximum effect \eqn{0.5 * E_{max}}.
##'
##' All parameters can currently be assigned uniform distributions or be fixed
##' in advance.
##'
##' @slot E0 either a fixed number or the two uniform distribution parameters
##' @slot Emax either a fixed number or the two uniform distribution parameters
##' @slot ED50 either a fixed number or the two uniform distribution parameters
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-DualEndpointEmax.R
##' @export
##' @keywords classes
.DualEndpointEmax <-
    setClass("DualEndpointEmax",
             representation(E0="numeric",
                            Emax="numeric",
                            ED50="numeric",
                            refDose="numeric"),
             prototype(E0=c(0, 100),
                       Emax=c(0, 500),
                       ED50=c(0,500),
                       refDose=1000,
                       useFixed=
                           list(sigma2W=TRUE,
                                rho=TRUE,
                                E0=FALSE,
                                Emax=FALSE,
                                ED50=FALSE)),
             contains="DualEndpoint",
             validity=
                 function(object){
                     o <- Validate()

                     ## check the prior parameters with variable content
                     for(parName in c("E0", "Emax", "ED50"))
                     {
                         ## if we use a fixed value for this parameter
                         if(object@useFixed[[parName]])
                         {
                             ## check range of value
                             o$check(slot(object, parName) > 0,
                                     paste(parName, "must be positive"))
                         } else {
                             ## use a Uniform(a, b) prior
                             o$check(all(slot(object, parName) >= 0) &&
                                         (diff(slot(object, parName)) > 0),
                                     paste(parName,
                                           "has not proper prior parameters"))
                         }
                     }

                     ## check the refDose
                     o$check(object@refDose > 0,
                             "refDose must be positive")

                     o$result()
                 })
validObject(.DualEndpointEmax())

##' Initialization function for the "DualEndpointEmax" class
##'
##' @param E0 see \code{\linkS4class{DualEndpointEmax}}
##' @param Emax see \code{\linkS4class{DualEndpointEmax}}
##' @param ED50 see \code{\linkS4class{DualEndpointEmax}}
##' @param refDose see \code{\linkS4class{DualEndpointEmax}}
##' @param \dots additional parameters, see \code{\linkS4class{DualEndpoint}}
##' @return the \code{\linkS4class{DualEndpointEmax}} object
##'
##' @export
##' @keywords methods
DualEndpointEmax <- function(E0,
                             Emax,
                             ED50,
                             refDose,
                             ...)
{
    ## call the initialize function from DualEndpoint
    ## to get started
    start <- DualEndpoint(...)

    ## we need the dose grid here in the BUGS model,
    ## therefore add it to datanames
    start@datanames <- c(start@datanames,
                         "doseGrid")

    ## Find out which of the additional parameters are fixed
    for(parName in c("E0", "Emax", "ED50"))
    {
        start@useFixed[[parName]] <-
            identical(length(get(parName)), 1L)
    }

    ## build together the prior model and the parameters
    ## to be saved during sampling
    ## ----------

    start@priormodel <-
        joinModels(start@priormodel,
                   function(){

                       for (j in 1:nGrid)
                       {
                           StandDose[j] <- doseGrid[j] / refDose
                           betaW[j] <- E0 + (Emax - E0) * StandDose[j] /
                                            (ED50 + StandDose[j])
                       }
                   })

    ## we will fill in more, depending on which parameters
    ## are fixed, in these two variables:
    start@sample <- c(start@sample,
                      "betaW")
    newInits <- list()
    newModelspecs <- list(refDose=refDose)

    ## for E0:
    if(! start@useFixed[["E0"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for E0
                           E0 ~ dunif(E0low, E0high)
                       })

        start@sample <- c(start@sample,
                          "E0")

        newInits$E0 <- mean(E0)
        newModelspecs$E0low <- E0[1]
        newModelspecs$E0high <- E0[2]
    } else {
        newModelspecs$E0 <- E0
    }

    ## for Emax:
    if(! start@useFixed[["Emax"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for Emax
                           Emax ~ dunif(EmaxLow, EmaxHigh)
                       })

        start@sample <- c(start@sample,
                          "Emax")

        newInits$Emax <- mean(Emax)
        newModelspecs$EmaxLow <- Emax[1]
        newModelspecs$EmaxHigh <- Emax[2]
    } else {
        newModelspecs$Emax <- Emax
    }

    ## for ED50:
    if(! start@useFixed[["ED50"]])
    {
        start@priormodel <-
            joinModels(start@priormodel,
                       function(){
                           ## uniform for ED50
                           ED50 ~ dunif(ED50Low, ED50High)
                       })

        start@sample <- c(start@sample,
                          "ED50")

        newInits$ED50 <- mean(ED50)
        newModelspecs$ED50Low <- ED50[1]
        newModelspecs$ED50High <- ED50[2]
    } else {
        newModelspecs$ED50 <- ED50
    }


    ## now define the new modelspecs and init functions:
    oldModelspecs <- start@modelspecs
    start@modelspecs <- function()
    {
        c(oldModelspecs(),
          newModelspecs)
    }

    oldInit <- start@init
    start@init <- function(y, w, nGrid)
    {
        c(oldInit(y, w, nGrid),
          newInits)
    }

    ## finally call the constructor
    .DualEndpointEmax(start,
                      E0=E0,
                      Emax=Emax,
                      ED50=ED50,
                      refDose=refDose)
}
validObject(DualEndpointEmax(E0=10,
                             Emax=50,
                             ED50=20,
                             refDose=10,
                             mu=c(0, 1),
                             Sigma=diag(2),
                             sigma2W=1,
                             rho=0))


## ============================================================



##' Standard logistic model with flexible mixture of two bivariate normal priors
##'
##' This is standard logistic regression model with a mixture of two bivariate
##' normal priors on the intercept and slope parameters. The weight of the two
##' normal priors is a model parameter, hence it is a flexible mixture.
##' This type of prior is often used with a mixture of a minimal informative
##' and an informative component, in order to make the CRM more robust to
##' data deviations from the informative component.
##'
##' The covariate is the natural logarithm of the dose \eqn{x} divided by
##' the reference dose \eqn{x^{*}}:
##'
##' \deqn{logit[p(x)] = \alpha + \beta \cdot \log(x/x^{*})}
##' where \eqn{p(x)} is the probability of observing a DLT for a given dose
##' \eqn{x}.
##'
##' The prior is
##' \deqn{(\alpha, \beta) \sim
##' w * Normal(\mu_{1}, \Sigma_{1}) + (1 - w) * Normal(\mu_{2}, \Sigma_{2})}
##'
##' The weight w for the first component is assigned a beta prior B(a, b).
##'
##' The slots of this class comprise two lists, containing the mean vector, the
##' covariance and precision matrices of the two bivariate normal distributions
##' each, the parameters of the beta prior for the first component weight, as
##' well as the reference dose.
##'
##' @slot comp1 the specifications of the first component: a list with
##' \code{mean}, \code{cov} and \code{prec} for the first bivariate normal prior
##' @slot comp2 the specifications of the second component
##' @slot weightpar the beta parameters for the weight of the first component
##' @slot refDose the reference dose \eqn{x^{*}}
##'
##' @example examples/Model-class-LogisticNormalMixture.R
##' @export
##' @keywords classes
.LogisticNormalMixture <-
    setClass(Class="LogisticNormalMixture",
             contains="Model",
             representation(comp1="list",
                            comp2="list",
                            weightpar="numeric",
                            refDose="numeric"),
             prototype(comp1=
                           list(mean=c(0, 1),
                                cov=diag(2),
                                prec=diag(2)),
                       comp2=
                           list(mean=c(-1, 1),
                                cov=diag(2),
                                prec=diag(2)),
                       weightpar=c(a=1, b=1),
                       refDose=1),
             validity=
                 function(object){
                     o <- Validate()

                     for(thisSlot in c("comp1", "comp2"))
                     {
                         thisList <- slot(object, thisSlot)
                         o$check(identical(names(thisList),
                                           c("mean", "cov", "prec")),
                                 paste(thisSlot,
                                       "must be a list with elements mean, cov, prec"))
                         o$check(length(thisList$mean) == 2,
                                 paste(thisSlot,
                                       "mean must have length 2"))
                         o$check(identical(dim(thisList$cov), c(2L, 2L)) &&
                                     ! is.null(chol(thisList$cov)),
                                 paste(thisSlot,
                                       "cov must be positive-definite 2x2 covariance matrix"))
                         o$check(all.equal(solve(thisList$cov), thisList$prec,
                                           check.attributes=FALSE),
                                 paste(thisSlot,
                                       "prec must be inverse of cov"))
                     }

                     o$check(all(object@weightpar > 0) &&
                                 identical(names(object@weightpar),
                                           c("a", "b")),
                             "weightpar does not specify proper prior parameters")
                     o$check(is.scalar(object@refDose) && (object@refDose > 0),
                             "refDose must be scalar and positive")
                 })
validObject(.LogisticNormalMixture())


##' Initialization function for the "LogisticNormalMixture" class
##'
##' @param comp1 the specifications of the first component: a list with
##' \code{mean} and \code{cov} for the first bivariate normal prior
##' @param comp2 the specifications of the second component
##' @param weightpar the beta parameters for the weight of the first component
##' @param refDose the reference dose
##' @return the \code{\linkS4class{LogisticNormalMixture}} object
##'
##' @export
##' @keywords methods
LogisticNormalMixture <- function(comp1,
                                  comp2,
                                  weightpar,
                                  refDose)
{
    ## add precision matrices to component lists
    comp1 <- c(comp1,
               list(prec=solve(comp1$cov)))
    comp2 <- c(comp2,
               list(prec=solve(comp2$cov)))

    .LogisticNormalMixture(comp1=comp1,
                           comp2=comp2,
                           weightpar=weightpar,
                           refDose=refDose,
                           datamodel=
                               function(){
                                   ## the logistic likelihood:
                                   ## not changed from non-mixture case
                                   for (i in 1:nObs)
                                   {
                                       y[i] ~ dbern(p[i])
                                       logit(p[i]) <- alpha0 + alpha1 * StandLogDose[i]
                                       StandLogDose[i] <- log(x[i] / refDose)
                                   }
                               },
                           priormodel=
                               function(){
                                   ## the multivariate normal prior on the coefficients
                                   theta[1:2] ~ dmnorm(priorMean[1:2, comp],
                                                       priorPrec[1:2, 1:2, comp])
                                   ## this is conditional on the component index
                                   ## "comp"

                                   ## component index is 1 or 2
                                   comp <- comp0 + 1

                                   ## it is 1 with probability w and
                                   ## 2 with probability 1 - w
                                   comp0 ~ dbern(wc)
                                   wc <- 1 - w

                                   ## we have a beta prior on w
                                   w ~ dbeta(weightpar[1], weightpar[2])

                                   ## extract actual coefficients
                                   alpha0 <- theta[1]
                                   alpha1 <- theta[2]

                                   ## dummy to use refDose here.
                                   ## It is contained in the modelspecs list below,
                                   ## so it must occur here
                                   bla <- refDose + 1
                               },
                           datanames=c("nObs", "y", "x"),
                           modelspecs=
                               function(){
                                   list(refDose=refDose,
                                        priorMean=
                                            cbind(comp1$mean,
                                                  comp2$mean),
                                        priorPrec=
                                            array(data=
                                                      c(comp1$prec,
                                                        comp2$prec),
                                                  dim=c(2, 2, 2)),
                                        weightpar=weightpar)
                               },
                           dose=
                               function(prob, alpha0, alpha1){
                                   StandLogDose <- (logit(prob) - alpha0) / alpha1
                                   return(exp(StandLogDose) * refDose)
                               },
                           prob=
                               function(dose, alpha0, alpha1){
                                   StandLogDose <- log(dose / refDose)
                                   return(plogis(alpha0 + alpha1 * StandLogDose))
                               },
                           init=
                               ## todo: find better starting values
                               function(){
                                   list(theta=c(0, 1))
                               },
                           sample=
                               c("alpha0", "alpha1", "w"))
}
validObject(LogisticNormalMixture(comp1=
                                      list(mean=c(0, 1),
                                           cov=diag(2)),
                                  comp2=
                                      list(mean=c(-1, 1),
                                           cov=diag(2)),
                                  weightpar=c(a=1, b=1),
                                  refDose=1))

## ============================================================

##' Standard logistic model with online mixture of two bivariate log normal priors
##'
##' This model can be used when data is arising online from the informative
##' component of the prior, at the same time with the data of the trial of
##' main interest. Formally, this is achieved by assuming that the probability
##' of a DLT at dose \eqn{x} is given by
##' 
##' \deqn{p(x) = \pi p_{1}(x) + (1 - \pi) p_{2}(x)}
##' 
##' where \eqn{\pi} is the probability for the model \eqn{p(x)} being the same
##' as the model \eqn{p_{1}(x)} - this is 
##' the informative component of the prior. From this model data arises in 
##' parallel: at doses \code{xshare}, DLT information \code{yshare} is observed, 
##' in total \code{nObsshare} data points, see \code{\linkS4class{DataMixture}}.
##' On the other hand, \eqn{1 - \pi}
##' is the probability of a separate model \eqn{p_{2}(x)}. Both components 
##' have the same log normal prior distribution, which can be specified by the
##' user, and which is inherited from the \code{\linkS4class{LogisticLogNormal}}
##' class.
##' 
##' @slot shareWeight the prior weight for sharing the same model \eqn{p_{1}(x)} 
##'
##' @seealso the \code{\linkS4class{DataMixture}} class for use with this model
##' @example examples/Model-class-LogisticLogNormalMixture.R
##' @export
##' @keywords classes
.LogisticLogNormalMixture <-
  setClass(Class="LogisticLogNormalMixture",
           contains="LogisticLogNormal",
           representation(shareWeight="numeric"),
           prototype(shareWeight=0.1),
           validity=
             function(object){
               o <- Validate()
               
               o$check(is.probability(object@shareWeight),
                       "shareWeight does not specify a probability")
               
             })
validObject(.LogisticLogNormalMixture())


##' Initialization function for the "LogisticLogNormalMixture" class
##'
##' @param mean the prior mean vector
##' @param cov the prior covariance matrix
##' @param refDose the reference dose
##' @param shareWeight the prior weight for the share component
##' @return the \code{\linkS4class{LogisticLogNormalMixture}} object
##'
##' @export
##' @keywords methods
LogisticLogNormalMixture <- function(mean,
                                     cov,
                                     refDose,
                                     shareWeight)
{
  .LogisticLogNormalMixture(mean=mean,
                            cov=cov,
                            refDose=refDose,
                            shareWeight=shareWeight,
                            datamodel=
                              function(){
                                ## the logistic likelihood:
                                
                                ## mixture for the new combo obs
                                for (i in 1:nObs)
                                {
                                  ## the bernoulli distribution:
                                  y[i] ~ dbern(p[comp, i])
                                  
                                  ## comp gives the component -
                                  ## non-informative (1) or share (2)
                                  
                                  ## the two components:
                                  for (k in 1:2)
                                  {
                                    logit(p[k, i]) <- alpha0[k] + alpha1[k] *
                                      StandLogDose[i]
                                  }
                                  
                                  ## just the standardized log dose:
                                  StandLogDose[i] <- log(x[i] / refDose)
                                }
                                
                                ## just from share for the share obs
                                for (j in 1:nObsshare)
                                {
                                  ## the bernoulli distribution:
                                  yshare[j] ~ dbern(pshare[j])
                                  
                                  ## take the correct - second - component
                                  logit(pshare[j]) <- alpha0[2] + alpha1[2] *
                                    StandLogDoseshare[j]
                                  
                                  ## just the standardized log dose:
                                  StandLogDoseshare[j] <- log(xshare[j] / refDose)
                                }
                                
                              },
                            priormodel=
                              function(){
                                
                                ## compute precision matrix
                                priorPrec[1:2,1:2] <- inverse(cov[,])
                                
                                ## the two components: same prior
                                for (k in 1:2)
                                {
                                  theta[k, 1:2] ~ dmnorm(mean[1:2],
                                                         priorPrec[1:2,1:2])
                                  
                                  alpha0[k] <- theta[k, 1]
                                  alpha1[k] <- exp(theta[k, 2])
                                }
                                
                                ## the component indicator
                                comp ~ dcat(catProbs)
                                
                                ## dummy to use refDose here.
                                ## It is contained in the modelspecs list below,
                                ## so it must occur here
                                bla <- refDose + 1
                              },
                            datanames=c("nObs", "y", "x", "nObsshare", "yshare", "xshare"),
                            modelspecs=
                              function(){
                                list(mean=mean,
                                     cov=cov,
                                     refDose=refDose,
                                     catProbs=c(1 - shareWeight, shareWeight))
                              },
                            dose=
                              function(prob, alpha0, alpha1, comp){
                                stop("not implemented")
                              },
                            prob=
                              function(dose, alpha0, alpha1, comp){
                                StandLogDose <- log(dose / refDose)
                                selectMat <- cbind(seq_len(nrow(alpha0)), comp)
                                return(plogis(alpha0[selectMat] +
                                                alpha1[selectMat] * StandLogDose))
                              },
                            init=
                              function(){
                                list(theta=matrix(c(0, 0, 1, 1), nrow=2))
                              },
                            sample=
                              c("alpha0", "alpha1", "comp"))
}
validObject(LogisticLogNormalMixture(mean=c(0, 1),
                                     cov=diag(2),
                                     shareWeight=0.1,
                                     refDose=1))


## ============================================================

##' Standard logistic model with fixed mixture of multiple bivariate (log) normal priors
##'
##' This is standard logistic regression model with a mixture of multiple bivariate
##' (log) normal priors on the intercept and slope parameters. The weights of the
##' normal priors are fixed, hence no additional model parameters are introduced.
##' This type of prior is often used to better approximate a given posterior
##' distribution, or when the information is given in terms of a mixture.
##'
##' The covariate is the natural logarithm of the dose \eqn{x} divided by
##' the reference dose \eqn{x^{*}}:
##'
##' \deqn{logit[p(x)] = \alpha + \beta \cdot \log(x/x^{*})}
##' where \eqn{p(x)} is the probability of observing a DLT for a given dose
##' \eqn{x}.
##'
##' The prior is
##' \deqn{(\alpha, \beta) \sim
##' \sum_{j=1}^{K} w_{j} Normal(\mu_{j}, \Sigma_{j})}
##' if a normal prior is used and
##' \deqn{(\alpha, \log(\beta)) \sim
##' \sum_{j=1}^{K} w_{j} Normal(\mu_{j}, \Sigma_{j})}
##' if a log normal prior is used.
##'
##' The weight \eqn{w_{j}} of the components are fixed and sum to 1.
##'
##' The (additional) slots of this class comprise two lists, containing the mean
##' vector, the covariance and precision matrices of the two bivariate normal
##' distributions each, the parameters of the beta prior for the first component
##' weight, as well as the reference dose. Moreover, a slot specifies whether a
##' log normal prior is used.
##'
##' @slot components a list with one entry per component of the mixture.
##' Each entry is a list with \code{mean}, \code{cov} and \code{prec} for the
##' bivariate normal prior
##' @slot weights the weights of the components, these must be positive and sum
##' to 1
##' @slot refDose the reference dose \eqn{x^{*}}
##' @slot logNormal is a log normal prior specified for each of the components?
##'
##' @example examples/Model-class-LogisticNormalFixedMixture.R
##' @export
##' @keywords classes
.LogisticNormalFixedMixture <-
    setClass(Class="LogisticNormalFixedMixture",
             contains="Model",
             representation(components="list",
                            weights="numeric",
                            refDose="numeric",
                            logNormal="logical"),
             prototype(components=
                           list(comp1=
                                    list(mean=c(0, 1),
                                         cov=diag(2),
                                         prec=diag(2)),
                                comp2=
                                    list(mean=c(-1, 1),
                                         cov=diag(2),
                                         prec=diag(2))),
                       weights=c(1/2, 1/2),
                       refDose=1,
                       logNormal=TRUE),
             validity=
                 function(object){
                     o <- Validate()

                     for(thisComp in seq_along(object@components))
                     {
                         thisList <- object@components[[thisComp]]
                         o$check(identical(names(thisList),
                                           c("mean", "cov", "prec")),
                                 paste("element", thisComp,
                                       "must be a list with elements mean, cov, prec"))
                         o$check(length(thisList$mean) == 2,
                                 paste("element", thisComp,
                                       "mean must have length 2"))
                         o$check(identical(dim(thisList$cov), c(2L, 2L)) &&
                                     ! is.null(chol(thisList$cov)),
                                 paste("element", thisComp,
                                       "cov must be positive-definite 2x2 covariance matrix"))
                         o$check(all.equal(solve(thisList$cov), thisList$prec,
                                           check.attributes=FALSE),
                                 paste("element", thisComp,
                                       "prec must be inverse of cov"))
                     }

                     o$check(identical(length(object@components),
                                       length(object@weights)),
                             "components must have same length as weights")
                     o$check(all(object@weights > 0) &&
                                 (sum(object@weights) == 1),
                             "weights must be positive and sum to 1")
                     o$check(is.scalar(object@refDose) && (object@refDose > 0),
                             "refDose must be scalar and positive")
                     o$check(is.bool(object@logNormal),
                             "logNormal must be TRUE or FALSE")
                 })
validObject(.LogisticNormalFixedMixture())

##' Initialization function for the "LogisticNormalFixedMixture" class
##'
##' @param components the specifications of the mixture components: a list with
##' one list of \code{mean} and \code{cov} for each bivariate (log) normal prior
##' @param weights the weights of the components, these must be positive and
##' will be normalized to sum to 1
##' @param refDose the reference dose
##' @param logNormal should a log normal prior be specified, such that the mean
##' vectors and covariance matrices are valid for the intercept and log slope?
##' (not default)
##' @return the \code{\linkS4class{LogisticNormalFixedMixture}} object
##'
##' @export
##' @keywords methods
LogisticNormalFixedMixture <- function(components,
                                       weights,
                                       refDose,
                                       logNormal=FALSE)
{
    ## add precision matrices to component lists
    components <- lapply(components,
                         function(x){
                             c(x,
                               list(prec=solve(x$cov)))})

    ## normalize the weights to sum to 1
    weights <- weights / sum(weights)

    ## go to the general initialize method now
    .LogisticNormalFixedMixture(components=components,
                                weights=weights,
                                refDose=refDose,
                                logNormal=logNormal,
                                datamodel=
                                    function(){
                                        ## the logistic likelihood:
                                        ## not changed from non-mixture case
                                        for (i in 1:nObs)
                                        {
                                            y[i] ~ dbern(p[i])
                                            logit(p[i]) <- alpha0 + alpha1 * StandLogDose[i]
                                            StandLogDose[i] <- log(x[i] / refDose)
                                        }
                                    },
                                priormodel=
                                    if(logNormal){
                                        function()
                                        {
                                            ## the multivariate normal prior on the coefficients
                                            theta[1:2] ~ dmnorm(priorMean[1:2, comp],
                                                                priorPrec[1:2, 1:2, comp])
                                            ## this is conditional on the component index
                                            ## "comp"

                                            ## mixture for component index
                                            comp ~ dcat(weights)

                                            ## extract actual coefficients
                                            alpha0 <- theta[1]
                                            alpha1 <- exp(theta[2])
                                            ## single difference to !logNormal ...

                                            ## dummy to use refDose here.
                                            ## It is contained in the modelspecs list below,
                                            ## so it must occur here
                                            bla <- refDose + 1
                                        }
                                    } else {
                                        function()
                                        {
                                            ## the multivariate normal prior on the coefficients
                                            theta[1:2] ~ dmnorm(priorMean[1:2, comp],
                                                                priorPrec[1:2, 1:2, comp])
                                            ## this is conditional on the component index
                                            ## "comp"

                                            ## mixture for component index
                                            comp ~ dcat(weights)

                                            ## extract actual coefficients
                                            alpha0 <- theta[1]
                                            alpha1 <- theta[2]

                                            ## dummy to use refDose here.
                                            ## It is contained in the modelspecs list below,
                                            ## so it must occur here
                                            bla <- refDose + 1
                                        }
                                    },
                                datanames=c("nObs", "y", "x"),
                                modelspecs=
                                    function(){
                                        list(refDose=refDose,
                                             priorMean=
                                                 do.call(cbind,
                                                         lapply(components, "[[", "mean")),
                                             priorPrec=
                                                 array(data=
                                                           do.call(c,
                                                                   lapply(components, "[[", "prec")),
                                                       dim=c(2, 2, length(components))),
                                             weights=weights)
                                    },
                                dose=
                                    function(prob, alpha0, alpha1){
                                        StandLogDose <- (logit(prob) - alpha0) / alpha1
                                        return(exp(StandLogDose) * refDose)
                                    },
                                prob=
                                    function(dose, alpha0, alpha1){
                                        StandLogDose <- log(dose / refDose)
                                        return(plogis(alpha0 + alpha1 * StandLogDose))
                                    },
                                init=
                                    ## todo: find better starting values
                                    function(){
                                        list(theta=c(0, 1))
                                    },
                                sample=
                                    c("alpha0", "alpha1"))
}
validObject(LogisticNormalFixedMixture(components=
                                           list(comp1=
                                                    list(mean=c(0, 1),
                                                         cov=diag(2),
                                                         prec=diag(2)),
                                                comp2=
                                                    list(mean=c(-1, 1),
                                                         cov=diag(2),
                                                         prec=diag(2))),
                                       weights=c(1/2, 1/2),
                                       refDose=1))

## =========================================================================
##' Class of models using expressing their prior in form of Pseudo data
##' 
##' This is the Pseudo model class, from which all models where their prior 
##' are expressed in form of pseudo data (as if some data are 
##' available before the trial starts) inherit. It also inherits all slots
##' from \code{\linkS4class{AllModels}}.No slots for this class
##' 
##' @seealso \code{\linkS4class{LogisticIndepBeta}},
##' \code{\linkS4class{Effloglog}},
##' \code{\linkS4class{EffFlexi}}
##'  
##' @export
##' @keywords classes
.ModelPseudo<-setClass(Class="ModelPseudo",
                       contains="AllModels"
)
validObject(.ModelPseudo)
##' No intialization function

## ===========================================================================

##' Class for DLE models using pseudo data prior. 
##' This is a class of DLE (dose-limiting events) models/ toxicity model which contains all DLE models 
##' for which their prior are specified in form of pseudo data (as if there is some data before
##' the trial starts). It inherits all slots from \code{\linkS4class{ModelPseudo}}
##' 
##' The \code{dose} function has a first argument \code{prob}, a scalar a probability of 
##' the occurrence of a DLE which is targeted. Additional arguments are models parameters. 
##' It computes, using the model parameter(s)/ model parameter(s) samples, the resulting dose. 
##' Note that the model parameters are called exactly as in the \code{model}. The model estimates 
##' generated can be single values of the maximum likelihodd estimates (prior or posterior modal
##' estimates) or samples of the model estimates generated. If samples of the model estimates are
##' generated, the model parameters (samples) must be included in the \code{samples} vector.
##' The vectors of all samples for these model paramters will be supplied to the function such 
##' that the function will be able to process vectors of model parameters.
##' 
##' The \code{prob} function has a first argument \code{dose}, a scalar dose level which is targeted.
##' Additional arguments are model paramters. It computes using model paramter(s) (samples), the 
##' resulting probabilities of a DLE occuring at the target dose level. If samples of model parameters
##' are generated, the function must vectorize over the model parameters.
##' 
##' Note that \code{dose} and \code{prob} are the inverse functions of each other.
##' 
##' The \code{data} must obey the covention that the data input is called exactly in the 
##' \code{\linkS4class{Data}} class. This refers to any observed DLE responses (\code{y} in 
##' \code{\linkS4class{Data}} class), the dose (levels) (\code{x} in \code{\linkS4class{Data}} class)
##' at which these responses are observed, all dose levels considered in the study (\code{doseGrid}
##' in \code{\linkS4class{Data}}) class and other specifications in \code{\linkS4class{Data}}
##' class that can be used to generate prior or
##' posterior modal estimates or samples estimates for model parmater(s). If no responses is observed,
##' at least \code{doseGrid} in \code{\linkS4class{Data}} has to be specified in \code{data} slot for which
##' prior modal estimates or samples can be obtained for model parameters based on the specified pseudo 
##' data.
##' 
##' 
##' 
##' @slot dose a function computing the dose level reaching a specific target probabilty of the occurrence 
##' of a DLE, based on the model parameters. The model paramters (samples)are obtained based on the prior 
##' specified in form of pseudo data and together with (if any) the observed
##' DLE responses and their corresponding dose levels (see details above)
##' @slot prob a function computing the probability of the occurrence of a DLEat a specidfied dose level, 
##' based on the model parameters. The model paramters (samples) are obtained the prior specified in form 
##' of pseudo data and together with (if any) the observed DLE responses and their 
##' corresponding dose levels (see dtails above)
##' @slot data refers to the data input specification in \code{\linkS4class{Data}} class which are used to
##' obtain model paramters estimates or samples (see details above)
##'
##' 
##' @seealso \code{\linkS4class{LogisticIndepBeta}},
##' \code{\linkS4class{Effloglog}},
##' \code{\linkS4class{EffFlexi}}
##' 
##' 
##' @export
##' @keywords classes
.ModelTox<-setClass(Class="ModelTox",
                    representation(dose="function",
                                   prob="function",
                                   data="Data"),
                    contains="ModelPseudo"
)
validObject(.ModelTox)
##' No Initialization function

## ==========================================================================================

##' class for Efficacy models using pseudo data prior
##' 
##' This is a class of which contains all efficacy models for which their prior are specified in 
##' form of pseudo data. It inherits all slots from \code{\linkS4class{ModelPseudo}}
##' 
##' The \code{dose} function has a first argument \code{ExpEff}, a scalar expected efficacy value 
##' which is targeted. Additional arguements are model parameters. It computes using modal estimate(s)
##' or samples model parameter(s), the resulting expected efficacy value at that dose level. If samples
##' of the model parameters are used, the function must vectorize over the model parameters.
##' 
##' The \code{ExpEff} function has a first argument \code{dose}, a scalar dose level which is targeted. 
##' Additional arguments are model parameters. It computes using modal estimates or samples of the
##' model parameter(s), the resulting dose level given that particular expected efficacy value. If samples
##' of the model parameter(s) are used, the function must vectorize over the model parameters.
##' 
##' The \code{data} must obey the covention that the data input is called exactly in the 
##' \code{\linkS4class{DataDual}} class. This refers to any observed Efficacy/biomarker responses 
##' (\code{w} in 
##' \code{\linkS4class{DataDual}} class), the dose (levels) (\code{x} in \code{\linkS4class{DataDual}} or
##' \code{Data} class)
##' at which these responses are observed, all dose levels considered in the study (\code{doseGrid}
##' in \code{\linkS4class{DataDual}} or \code{Data}) class and other specifications in 
##' \code{\linkS4class{DataDual}}
##' class that can be used to generate prior or
##' posterior modal estimates or samples estimates for model parmater(s). If no responses is observed,
##' at least \code{doseGrid} in \code{\linkS4class{DataDual}} has to be specified in \code{data} slot
##' for which prior modal estimates or samples can be obtained for model parameters based on 
##' the specified pseudo data.
##' 
##' @slot dose a function computing the dose reaching a specific target value of expected efficacy, based
##' on the model parameter(s). The model parameter(s) (samples) are obtained based on prior specified 
##' in form of pseudo data and if any together with any observed responses (see details above)
##' 
##' @slot ExpEff a function computing the expected efficacy (value) for a specific dose, based on model 
##' parameter(s). The model parameter(s) (samples) are obtained based on pseudo data prior and (if any) 
##' with observed responses (see details above)
##' 
##' @slot data refers to the data input specification in \code{\linkS4class{DataDual}} class which are used to
##' obtain model paramters estimates or samples (see details above)
##'
##' @seealso \code{\linkS4class{Effloglog}},
##' \code{\linkS4class{EffFlexi}}
##' 
##' @export
##' @keywords classes
.ModelEff<-setClass(Class="ModelEff",
                    representation(dose="function",
                                   ExpEff="function",
                                   data="DataDual"),
                    contains="ModelPseudo"
)
validObject(.ModelEff)
##' No initialization function

## ==============================================================================
##' Standard logistic model with prior in form of pseudo data
##' 
##' This is a class for the two-parameter logistic regression DLE model with prior expressed
##' in form of pseudo data. This model describe the relationhship of the binary DLE (dose-limiting
##' events) responses and the dose levels. More specifically, this DLE model reprsents the relationship 
##' of the probabilities of the occurrence of a DLE with their corresponding dose levels in log scale.
##' This model is specified as 
##' \deqn{p(d_{(j)})= \frac{exp(\phi_1+\phi_2 log(d_{(j)}))}{1+exp(\phi_1+\phi_2 log(d_{(j)}))}}
##' for any dose j where \eqn{p(d_{(j)})} is the probability of the occurrence of a DLE at dose j.
##' The two parameters of this model is the intercept \eqn{\phi_1} and the slope \eqn{\phi_2}
##' It inherits all slots from \code{\linkS4class{ModelTox}} class.
##' 
##' The pseudo data can be interpreted as as if we obtain some observations before the trial starts.
##' These pseudo data can be used to express our prior, the initial beliefs for the model parameter(s).
##' The pseudo data are expressed in the following way. First, fix at least two dose levels which are
##' Then ask for experts' opinion how many subjects are to be treated at each of these dose levels and
##' the number of subjects observed with DLE are observed. At each dose level, the number of subjects 
##' observed with a DLE divided by the total number of subjects treated is the probability of the 
##' occurrence of a DLE at that particular dose level. The probabilities of the occurrence of a DLE
##' based on these pseudo data are independent Beta distributions. Therefore, the joint prior probability 
##' density function of all these probabilities can be obtained. Hence, by a change of variable, the 
##' joint prior probability density function of the two parameters in this model can also be obtained.
##' In addition, a conjugate joint perior density function of the two paramaters in the model is used.
##' For details about the form of all these joint prior and posterior probability density function, please 
##' refers to Whitehead and Willamson (1998). 
##' 
##' 
##' When expressing the pseudo data, \code{binDLE},\code{DLEdose} and \code{DLEweights} are used.
##' The \code{binDLE} represents the number of subjects observed with DLE. Note that, since the imaginary 
##' nature of the pseudo data, the number of subjects observed wtih DLE is not necesssary to be integer(s)
##' but any scalar value.
##' The \code{DLEdose} represents the dose levels at which the pseudo DLE responses (\code{binDLE}) are 
##' observed.
##' The \code{DLEweights} represents the total number of subjects treated.
##' Since at least two DLE pseudo responses are needed to obtain prior modal estimates (same as the maximum 
##' likelihood estimates) for the model parameters. \code{binDLE}, \code{DLEdose} and \code{DLEweights} must
##' all be vectors of at least length 2. Since given one pseudo DLE responses, the number of subjects observed 
##' with a DLE relates to at which dose level they are treated and the total number of of subjects treated at
##' this dose level. Therefore, each of the elements in any of the vectors of \code{binDLE}, \code{DLEdose} and 
##' \code{DLEweights} must have a corresponding elements in the other two vectors. A set of three values with
##' one of each in the vectors of \code{binDLE}, \code{DLEdose} and \code{DLEweights}. In this model, each of 
##' these three values must be specified in the same position as in each of the vector of \code{binDLE}, 
##' \code{DLEdose} and \code{DLEweights}. The order of the values or elements in one of the vector \code{binDLE}, 
##' \code{DLEdose} and \code{DLEweights} must corresponds to the values or elements specified in the other two 
##' vectors.
##' 
##' @slot binDLE represents the vector of pseudo DLE responses. This must be at least f length 2 and the 
##' order of its elements must corresponds to values specified in \code{DLEdose} and \code{DLEweights}. 
##' (see details from above)
##' @slot DLEdose represents the vector of the corresponding dose levels observed at each of the 
##' pseudo DLE responses (\code{binDLE}). This mus be at least of length 2 and the order of its elements
##' must corresponds to values specified in \code{binDLE} and \code{DLEweights}.
##' (see details from above)
##' @slot DLEweights refers to the total number of subjects treated at each of the pseudo dose level 
##' (\code{DLEdose}). This must be of length of at least 2 and the oreder of its elements must corresponds
##' to values specified in \code{binDLE} and \code{DLEdose}. (see details from above)
##' @slot phi1 refers the intercept of the model. This slot is used in output to display the resulting prior 
##' or posterior modal estimate of the intercept obtained based on the pseudo data and (if any) 
##' observed data/responses.
##' @slot phi2 refers to slope of the model. This slot is used in output to display the resulting prior or 
##' posterior modal estimate of the slope obtained based on the pseudo data and (if any) the observed data/responses.
##' @slot Pcov refers to the covariance matrix of the intercept (phi1) and the slope parameters (phi2) of the 
##' model. This is used in output to display the resulting prior and posterior covariance matrix of phi1 and 
##' phi2 obtained, based on the pseudo data and (if any) the observed data and responses. This slot is needed for 
##' internal purposes.
##'  
##' @example examples/Model-class-LogisticIndepBeta.R
##' @export
##' @keywords classes
.LogisticIndepBeta<-
  setClass(Class="LogisticIndepBeta",
           representation(binDLE="numeric",
                          DLEdose="numeric",
                          DLEweights="numeric",
                          
                          phi1="numeric",
                          phi2="numeric",
                          Pcov="matrix"),
           prototype(binDLE=c(0,0),
                     DLEdose=c(1,1),
                     DLEweights=c(1,1)),
           contains="ModelTox",
           validity=
             function(object){
               o <- Validate()
               ##Check if at least two pseudo DLE responses are given
               o$check(length(object@binDLE) >= 2,
                       "length of binDLE must be at least 2")
               ##Check if at least two weights for pseudo DLE are given 
               o$check(length(object@DLEweights) >= 2,
                       "length of DLEweights must be at least 2")
               ##Check if at least two corresponding dose levels are given for the pseudo DLE responses
               o$check(length(object@DLEdose) >= 2,
                       "length of DLEdose must be at least 2")
               ##Check if pseudo DLE responses have same length with it corresponding dose levels and weights
               o$check((length(object@binDLE)==length(object@DLEweights))&(length(object@binDLE)==length(object@DLEdose))&(length(object@DLEweights)==length(object@DLEdose)),
                       "length of binDLE, DLEweights, DLEDose must be equal")
               o$result()
             })
validObject(.LogisticIndepBeta())  

##' Intialization function for "LogisticIndepBeta" class
##' @param binDLE the number of subjects observed with a DLE, the pseudo DLE responses
##' @param DLEdose the corresponding dose levels for the pseudo DLE responses, pseudo dose levels
##' @param DLEweights the total number of subjects treated at each of the dose levels, pseudo weights
##' @param data the input data to update estimates of model parameters and 
##' follow the \code{\linkS4class{Data}} object class specification
##' @return the \code{\linkS4class{LogisticIndepBeta}}
##' 
##' @export
##' @keywords methods
LogisticIndepBeta <- function(binDLE,
                              DLEdose,
                              DLEweights,
                              data)
{##if no observed DLE(data)
  if (length(data@y)==0){
    w1<-DLEweights
    y1<-binDLE
    x1<-DLEdose} else {w1<-c(DLEweights,rep(1,data@nObs))
    ##combine pseudo and observed
    y1<-c(binDLE,data@y)
    x1<-c(DLEdose,data@x)}
  ##Fit the pseudo data and DLE responses with their corresponding dose levels
  FitDLE<-suppressWarnings(glm(y1/w1~log(x1),family=binomial(link="logit"),weights=w1))
  SFitDLE<-summary(FitDLE)
  ##Obtain parameter estimates for dose-DLE curve
  phi1<-coef(SFitDLE)[1,1]
  phi2<-coef(SFitDLE)[2,1]
  ## covariance matrix of phi1 and phi2
  Pcov <- vcov(FitDLE)
  
  .LogisticIndepBeta(binDLE=binDLE,
                     DLEdose=DLEdose,
                     DLEweights=DLEweights,
                     phi1=phi1,
                     phi2=phi2,
                     Pcov=Pcov,
                     datanames=c("nObs","y","x"),
                     data=data,
                     dose=function(prob,phi1,phi2){
                       LogDose<-((log(prob/(1-prob)))-phi1)/phi2
                       return(exp(LogDose))
                     },
                     prob=function(dose,phi1,phi2){
                       LogDose<-log(dose)
                       pj<-(exp(phi1+phi2*LogDose))/(1+exp(phi1+phi2*LogDose))
                       return(pj)
                     }
  )
}

## ======================================================================================================
##' Class for the linear log-log efficacy model using pseudo data prior
##' 
##' This is the efficacy model which describe the relationship of the continuous efficacy responses and 
##' the dose levels. More specifically, this is a model to describe the linear relationship between the 
##' continuous efficacy responses and its coressponding dose level in log-log scale. 
##' The efficacy log-log model is given as 
##' \deqn{y_i=\theta_1 +theta_2 log(log(d_i))+\epsilon_i}
##' where \eqn{y_i} is the efficacy responses
##' for subject i, \eqn{d_i} is the dose level treated for subject i and \eqn{\epsilon_i} is the random error 
##' term of efficacy model at subject i such that \eqn{\epsilon_i} has a normal distribution of mean 0 and 
##' variance \eqn{\sigma^2=\nu^{-1}}. This variance is assumed to be the same for all subjects.
##' 
##' There are three parameters in this model which is to intercept \eqn{\theta_1}, the slope \eqn{\theta_2} 
##' and the precision \eqn{\nu} of the efficay responses.
##' It inherit all slots from \code{\linkS4class{ModelEff}}
##' 
##' The prior of this model is specified in form of pseudo data. First at least two dose levels are fixed.
##' Then ask for experts' opinion about the efficacy values that can be obtained at each of the dose levels
##' if one subject is treated at each of these dose levels. The prior modal estimates (same as the maximum 
##' likelihood estimates) can be obtained for the intercept and slope paramters in this model.
##' 
##' The \code{Eff} and \code{Effdose} are used to represent the prior in form of the pseudo data. 
##' The \code{Eff} represents the pseudo scalar efficacy values. The \code{Effdose} represents the dose levels
##' at which these pseudo efficacy values are observed. These pseudo efficay values are always specified by
##' assuming one subject are treated in each of the dose levels. Since at least 2 pseudo efficacy values are 
##' needed to obtain modal estimates of the intercept and slope parameters, both \code{Eff} and \code{Effdose}
##' must be vector of at least length 2. The position of the values or elements specified in \code{Eff} or
##' \code{Effdose} must be corresponds to the same elements or values in the other vector.
##' 
##' The \code{nu} represents the prior presion \eqn{\nu} of the pseudo efficacy responses. It is also known as the inverse 
##' of the variance of the pseduo efficacy responses. The precision can be a fixed constant or having a gamma
##' distribution. Therefore, single scalar value, a fixed 
##' value of the precision can be specified. If not, two positive scalar values must be specified as the 
##' shape and rate parameter of the gamma distribution. If there are some observed efficacy responses available,
##' in the output, \code{nu} will display the updated value of the precision or the updated values for the 
##' parameters of the gamma distribution.
##' 
##' 
##' Given the variance of the pseudo efficacy responses, the joint prior distribution of the intercept \eqn{\theta_1}
##' (theta1) and the slope \eqn{\theta_2} (theta2) of this model is a bivariate normal distribution. 
##' A conjugate posterior joint distribution is also used for theta1 and theta2. The joint prior bivariate 
##' normal distribution has 
##' mean \eqn{\boldsymbol\mu_0} and covariance matrix \eqn{(\nu \mathbf{Q}_0)^{-1}}. \eqn{\boldsymbol\mu_0} is a 
##' \eqn{2 \times 1}
##' column vector contains the prior modal estimates of the intercept (theta1) and the slope (theta2). Based on 
##' \eqn{r} for \eqn{r \geq 2} pseudo efficacy responses specified, \eqn{\mathbf{X}_0} will be the 
##'\eqn{r \times 2} design matrix 
##' obtained for these pseudo efficacy responses. the matrix \eqn{\mathbf{Q}_0} will be calculated by 
##' \eqn{\mathbf{Q}_0=\mathbf{X}_0 \mathbf{X}^T_0} where \eqn{\nu} is the precision of the pseudo efficacy responses.
##' For the joint posterior bivariate distribution, we have \eqn{\boldsymbol{\mu}} as the mean and 
##' \eqn{(\nu\mathbf{Q}_0)^{-1}} as the covariance matrix. Here, \eqn{\boldsymbol\mu} is the column vector containing the 
##' posterior modal estimates
##' of the intercept (theta1) and the slope (theta2). The design matrix \eqn{\mathbf{X}} obtained based only on 
##' observed efficacy responses will give \eqn{\mathbf{Q}=\mathbf{X}\mathbf{X}^T} with \eqn{\nu} as the precision of 
##' the observed efficay responses. If no observed efficay responses are availble (i.e only pseudo 
##' efficay responses are used), the \code{vecmu}, \code{matX}, \code{matQ} and \code{vecY} represents 
##' \eqn{\boldsymbol\mu_0}, \eqn{\mathbf{X}_0}, \eqn{\mathbf{Q}_0} and the column vector of pseudo efficay responses,
##' respectively. If there are some observed efficacy responses, \code{vecmu}, \code{matX}, \code{matQ} 
##' and \code{vecY} will represent \eqn{\boldsymbol\mu}, \eqn{\mathbf{X}}, \eqn{\mathbf{Q}} and the column vector contains
##' all observed efficacy responses, respectively. (see details in about the form of prior and posterior distribution)
##' 
##' @slot Eff the pseudo efficacy response, the scalar efficacy values. This must be a vector of at least 
##' length 2. Each element or value here must represents responses treated based on one subject. The order
##'  of its elements must corresponds to the values presented in vector \code{Effdose} (see details above)
##' @slot Effdose the pseudo efficacy dose level. This is the dose levels at which the pseudo efficacy 
##' responses are observed at. This must be a vector of at least length 2 and the orde of its elements must
##' corresponds to values presented in vector \code{Eff} (see detial above)
##' @slot nu refers to the prior precision of pseudo efficacy responses. This is either a fixed value or a 
##' vector of elements \code{a}, a positive scalar for the shape parameter, and \code{b}, a positive scalar 
##' for the rate parameter for the gamma dsitribution. (see detail from above) 
##' @slot useFixed a logical value if \code{nu} specified is a fixed value or not. This slot is needed for 
##' internal purposes and not to be touched by the user.
##' @slot theta1 The intercept \eqn{\theta_1} parameter of this efficacy log-log model. This slot is used in output to display
##' the resulting prior or posterior modal estimates obtained based on the pseudo data and (if any) the 
##' observed data/ responses.
##' @slot theta2 The slope \eqn{theta_2} parameter of the efficacy log-lgo model. This slot is used in output to display 
##' the resulting prior or posterior modal estimates obtained based on the pseudo data and (if any) the 
##' observed data/ responses.
##' @slot Pcov refers to the covariance matrix of the intercept (phi1) and slope (phi2) paramters of this model.
##' This slot is used in output to display the covariance matrix obtained based on the pseudo data and (if any)
##' the observed data/responses. This slot is needed for internal purposes.
##' @slot vecmu is the column vector of the prior or the posterior modal estimates of the intercept (phi1) and 
##' the slope (phi2).
##' This slot is used in output to display as the mean of the prior or posterior bivariate normal distribtuion
##' for phi1 and phi2. (see details from above)
##' @slot matX is the design matrix based on either the pseudo or all observed efficacy response. This is used in 
##' output to display the design matrix for the pseudo or the observed efficacy responses (see details from above)
##' @slot matQ is the square matrix of multiplying the the design matrix with its transponse. This is represented 
##' either using the only the pseudo efficay responses or only with the observed efficacy responses. This is display 
##' in the output (see details from above)
##' @slot vecY is the column vector either contains the pseudo efficay responses or all the observed efficacy 
##' responses. This is used in output to display the pseudo or observed efficacy responses (see detail from above)
##' 
##' 
##'@example examples/Model-class-Effloglog.R
##'@export
##'@keywords methods
.Effloglog<-
  setClass(Class="Effloglog", 
           representation(Eff="numeric",
                          Effdose="numeric",
                          nu="numeric",
                          useFixed="logical",
                          theta1="numeric",
                          theta2="numeric",
                          Pcov="matrix",
                          vecmu="matrix",
                          matX="matrix",
                          matQ="matrix",
                          vecY="matrix"),
           prototype(Eff=c(0,0),
                     Effdose=c(1,1),
                     nu=1/0.025,
                     useFixed=TRUE),
           contains="ModelEff",
           validity=
             function(object){
               o <- Validate()
               
               o$check(length(object@Eff) >= 2,
                       "length of Eff must be at least 2")
               o$check(length(object@Effdose) >= 2,
                       "length of Effdose must be at least 2")
               o$check(length(object@Eff)==length(object@Effdose),
                       "length of Eff and Effdose must be equal")
               if (object@useFixed == "TRUE"){
                 o$check((length(object@nu)==1)&&(object@nu > 0),
                         "nu must be a single postive real number")} else {
                           o$check(identical(names(slot(object,"nu")),c("a","b")),
                                   "nu must have names 'a' and 'b' ")
                           o$check(all(slot(object,"nu") > 0),
                                   "nu must have positive prior paramters")
                           o$check(identical(length(object@nu),2L),
                                   "nu must have length at most 2")
                         }
               o$result()
             })
validObject(.Effloglog())

##' Initialization function for the "Effloglog" class
##' 
##' @param Eff the pseudo efficacy responses
##' @param Effdose the corresponding dose levels for the pseudo efficacy responses
##' @param nu the precision (inverse of the variance) of the efficacy responses
##' @param data the input data of \code{\linkS4class{DataDual}} class to update model estimates
##' @return the \code{\linkS4class{Effloglog}} object
##' 
##' @importFrom MASS ginv
##' @export
##' @keywords methods
Effloglog<-function(Eff,
                    Effdose,
                    nu,
                    data)
{if (!all(data@doseGrid >=1))
  stop("doseGrid in data must be greater or equal to 1 for Effloglog model")
  ##No observed Efficacy response
  if (length(data@w)==0){
    w1<-Eff
    x1<-Effdose} else {##Combine pseudo data and Observed Efficacy without DLE
      w1<-c(Eff,getEff(data)$wNoDLE)
      x1<-c(Effdose,getEff(data)$xNoDLE)
      w2<-getEff(data)$wNoDLE
      x2<-getEff(data)$xNoDLE} 
  
  
  ##Check if sigma2/nu is a fixed contant
  
  useFixed <- identical(length(nu), 1L)
  ##Fit pseudo and observed efficacy
  FitEff<-suppressWarnings(glm(w1~log(log(x1)),family=gaussian))
  SFitEff<-summary(FitEff)
  ##Obtain paramter estimates
  theta1<-coef(SFitEff)[1,1]
  theta2<-coef(SFitEff)[2,1]
  ##covariance matrix of theta1 and theta2
  Pcov <- vcov(FitEff)
  ##if sigma2/nu is not a fixed constant
  if (length(nu)==2){
    mu0<-matrix(c(theta1,theta2),2,1)
    vecmu<-mu0
    X0<-matrix(c(1,1,log(log(Effdose[1])),log(log(Effdose[2]))),2,2)
    matX<-X0
    Q0=t(X0)%*%X0
    matQ<-Q0
    vecY<-matrix(Eff,2,1)
    ##if there are some observed efficacy
    if (length(data@w)!=0){
      X<-matrix(c(rep(1,length(x2)), log(log(x2))), length(x2),2)
      matX<-X
      mu<-MASS::ginv(Q0+t(X)%*%X)%*%(Q0%*%mu0+t(X)%*%t(t(w2)))
      vecmu<-mu
      Q<-Q0+t(X)%*%X
      matQ<-Q
      vecY<-matrix(w2,length(w2),1)
      a<-nu[1]+(length(w2))/2
      b<-nu[2]+(t(w2)%*%t(t(w2))+t(mu0)%*%Q0%*%mu0-t(mu)%*%Q%*%mu)/2
      nu[1]<-a
      nu[2]<-b}} else {nu<-nu} 
  
  .Effloglog(Eff=Eff,
             Effdose=Effdose,
             nu=nu,
             useFixed=useFixed,
             datanames=c("nObs","w","x"),
             data=data,
             dose=function(ExpEff,theta1,theta2){
               LogDose<-exp((ExpEff-theta1)/theta2)
               return(exp(LogDose))},
             
             ExpEff=function(dose,theta1,theta2){
               
               return(theta1+theta2*log(log(dose)))},
             theta1=theta1,
             theta2=theta2,
             Pcov=Pcov,
             vecmu=vecmu,
             matX=matX,
             matQ=matQ,
             vecY=vecY
  )}

## =========================================================================================
##' Class for the efficacy model in flexible form for prior expressed in form of pseudo data
##' 
##' This is a class where a flexible form is used to describe the realtionship between the efficacy
##' responses and the dose levels. This flexible form aims to capture different shape for the 
##' dose-efficacy curve and the mean efficacy responses at each dose level are estimated using MCMC.
##' In addition, the first (RW1) or second order (RW2) random walk model can be used for smoothing data. That is 
##' the random walk model is used to model the first or the second order differnece of the mean 
##' efficacy responses to its neighbouring dose levels of their mean efficacy responses.
##' The flexible form is specified as 
##' \deqn{\mathbf{W}\vert\boldsymbol{\beta_w}, \sigma^2 \sim Normal (\mathbf{X}_w \boldsymbol{\beta_w}, \sigma^2 \mathbf{I})}
##' where \eqn{\mathbf{W}} represent the column vector of the efficacy responses, \eqn{\boldsymbol{\beta_w}}
##' is th column vector of the mean efficacy responses for all dose levels, \eqn{\mathbf{X_w}} is the 
##' design matrix with entries \eqn{I_{i(j)}} which gives a value 1 if subject i is allocated to 
##' dose j. The \eqn{\sigma^2} (sigma2) is the variance of the efficacy responses which can be either fixed or from
##' an inverse gamma distribution.
##' 
##' The RW1 model is given as 
##' \deqn{\beta_{W,(j)} - \beta_{W,(j-1)} \sim Normal(0, \sigma^{2}_{\beta_{W}})}
##' where \eqn{\beta_{W,(j)}} is the mean efficacy responses at dose j
##' For the RW2 is given as 
##' \deqn{\beta_{W,(j-2)} - 2 \beta_{W,(j-1)} + \beta_{W,(j)} \sim Normal(0, \sigma^{2}_{\beta_{W}})}
##' The variance parameter \eqn{\sigma^{2}_{\beta_{W}}}. The variance \eqn{\sigma^{2}_{\beta_{W}}} 
##' (sigma2betaW) will be the same at all dose levels and can
##' either be fixed or assigned an inverse gamma prior distribution.
##' 
##' The \code{Eff} and \code{Effdose} are the pseduo efficacy responses and dose levels at which these 
##' pseudo efficacy responses are observed at. (see more details for \code{\linkS4class{Effloglog}} class)
##' \code{Eff} and \code{Effdose} must be vector of at least length 2. The values or elements in vectors
##' \code{Eff} or \code{Effdose} must put in the same position with its corresponding value in the other 
##' vector. The \code{sigma2} is the prior variance of the flexible efficacy form. The variance is either specified
##' with a single scalar value (fixed) or postive scalar value have to be specified for the \code{a} shape and 
##' \code{b} slope parameter for th inverse gamme distribtuion. Similarly, \code{sigma2betaW} is the prior variance 
##' of the random walk model which can be specified with a single scalar (fixed) value or specifying positive 
##' scalar values for the shape \code{a} and rate \code{b} parameters for the inverse gamma distributions. 
##' This model will output the updated value or the updated values of the paramters of the inverse gamma
##' distributions for \eqn{sigma^2} (sigma2) and \eqn{\sigma^2_{\beta_W}} (sigma2betaW)
##' 
##' @slot Eff the pseudo efficacy responses. A vector of at least length 2 with the elements here and its 
##' corresponding value in \code{Effdose} must be specified in the same position. (see dtails above)
##' @slot Effdose the dose levels at which the pseudo efficacy responses are observed. This is a vector of at 
##' least length 2 and the elements here and its corrresponding value in \code{Eff} must be specified in the
##' same postion. (see details from above)
##' @slot sigma2 the prior variance of the flexible efficacy form. It can be specified with a single positive
##' scalar or specifying \code{a}, the shape and \code{b}, the rate parameter of the inverse gamma 
##' distribution. (see details from above)
##' @slot sigma2betaW the prior variance of the random walk model for the mean efficact responses. A single
##' positve scalar can be specified or specifying \code{a}, the shape and \code{b}, the rate parameter of 
##' the inverse gamma distribution (see details from above)
##' @slot useFixed a list of with logical value to each of the parameters \code{sigma2} and \code{sigma2betaw}
##' indicating whether a fixed value is used or not; this slot is needed for internal purposes and not to
##' be touched by the user.
##' @slot useRW1 for specifying the random walk model for the mean efficacy responses; if \code{TRUE}, 
##' first order random walk model is used, otherwise the second-order random walk model.
##' @slot designW is the design matrix for the efficacy responses. If only the pseudo efficacy responses 
##' are used, this will be the design matrix of the pseudo efficacy responses. If there are some observed
##' efficacy responses available. It will be the design matrix based on both the pseudo and the observed
##' efficacy responses. 
##' @slot RWmat is the the difference matrix for the random walk model. This slot is needed for internal 
##' purposes and not to be touched by the user.
##' @slot RWmatRank is the rank of the difference matrix. This slot is needed for internal purposes and not
##' to be touched by the user.
##'
##' @example examples/Model-class-EffFlexi.R
##' @export
##' @keywords class
.EffFlexi<-setClass(Class="EffFlexi",
representation(Eff="numeric",
               Effdose="numeric",
               sigma2="numeric",
               sigma2betaW="numeric",
               useFixed="list",
               useRW1="logical",
               designW="matrix",
               RWmat="matrix",
               RWmatRank="integer"),
prototype(Eff=c(0,0),
          Effdose=c(1,1),
          sigma2=0.025,
          sigma2betaW=1,
          useRW1=TRUE,
          useFixed=list(sigma2=TRUE,sigma2betaW=TRUE)),
contains="ModelEff",
validity=
  function(object){
    o<- Validate()
    o$check(length(object@Eff) >= 2,
            "length of Eff must be at least 2")
    o$check(length(object@Effdose) >= 2,
            "length of Effdose must be at least 2")
    o$check(length(object@Eff)==length(object@Effdose),
            "length of Eff and Effdose must be equal")
    for (parName in c("sigma2","sigma2betaW"))
    {
      if (object@useFixed[[parName]]){
        o$check(slot(object,parName) > 0,
                paste(parName, "must be positive"))} else {
                  o$check(identical(names(slot(object,parName)),c("a","b")), 
                          paste(parName,"must have names 'a' and 'b'"))
                  o$check(all(slot(object,parName) > 0),
                          paste(parName, "must have positive prior parameters"))
                }
    }
    o$result()
  })
validObject(.EffFlexi())

##' Initialization function for the "EffFlexi" class
##' 
##' @param Eff the pseudo efficacy responses
##' @param Effdose the corresponding dose levels for the pseudo efficacy responses
##' @param sigma2 the prior variance of the efficacy responses which can be specified 
##' with a single positive scalar or with two positive scalar values fot the shape \code{a} and 
##' the rate \code{b} parameters of the inverse gamma distribtuion.
##' @param sigma2betaW the prior variance of the random walk model used for smoothing which can be
##' specified with a single postive scalar or with two positive scalars representing the shape \code{a}
##' and the rate \code{b} parameter of the inverse gamma distribution.
##' @param smooth used for smoothing data for this efficacy model. That is either the "RW1", the 
##' first-order random walk model or "RW2", the second-order random walk model is used of the mean
##' efficacy responses.
##' @param data the input data to update estimates of model parameters and 
##' follow the \code{\linkS4class{DataDual}} object class specification
##' @return the \code{\linkS4class{EffFlexi}} class object
##' 
##' @export
##' @keywords methods

EffFlexi <- function(Eff,
                     Effdose,
                     sigma2,
                     sigma2betaW,
                     smooth=c("RW1","RW2"),
                     data
)
{##No observed Efficacy response
  if (length(data@w)==0){
    w1<-Eff
    x1<-Effdose} else {
      ## with observed efficacy responses and no DLE observed
      w1<-c(Eff,getEff(data)$wNoDLE)
      x1<-c(Effdose,getEff(data)$xNoDLE)
    }
  ## Match dose levels in x1 with the all dose levels for evaluations
  x1Level <- match(x1,data@doseGrid)
  smooth<-match.arg(smooth)
  useRW1<- smooth == "RW1"
  useFixed<-list()
  for (parName in c("sigma2","sigma2betaW"))
  {useFixed[[parName]] <- identical(length(get(parName)),1L)}
  #design matrics
  designW <- model.matrix(~ -1 + I(factor(x1Level, levels=seq_len(data@nGrid))))
  dimnames(designW) <- list(NULL,NULL)
  
  ##difference matrix of order 1:
  D1mat<- cbind(0,diag(data@nGrid-1)) - cbind(diag(data@nGrid - 1),0)
  
  ## set up the random walk penalty matrix and its rank:
  if (useRW1)
  {## the rank-deficient prior precision for the RW1 prior:
    RWmat <- crossprod(D1mat)
    ##Rank: dimension -1
    RWmatRank <- data@nGrid-1L
  } else {##second-order difference
    D2mat <- D1mat[-1,-1] %*% D1mat
    RWmat <- crossprod(D2mat)
    RWmatRank <- data@nGrid-2L
  }
  .EffFlexi(Eff=Eff,
            Effdose=Effdose,
            sigma2=sigma2,
            sigma2betaW=sigma2betaW,
            datanames=c("nObs","w","x"),
            data=data,
            dose=function(ExpEff){
              ##Find dose level given a particular Expected Efficacy level with linear Interpolation
              dosevec<-c()
              for (k in 1:sampleSize(options)){
                IxEff0<- max(which((ExpEff-Effsamples@data$ExpEff[k,]) >= 0))
                IxEff1<- min(which((ExpEff-Effsamples@data$ExpEff[k,]) < 0))
                Interpoldose<-data@doseGrid[IxEff0]+(data@doseGrid[IxEff1]-data@doseGrid[IxEff0])*((ExpEff-Effsamples@data$ExpEff[k,IxEff0])/(Effsamples@data$ExpEff[k,IxEff1]-Effsamples@data$ExpEff[k,IxEff0]))
                dosevec[k]<-Interpoldose
              }
              ##return coreresponding dose levels
              return(dosevec)},
            
            ExpEff=function(dose,data,Effsamples){
              ##Find the ExpEff with a given dose level
              ##Check if given dose is in doseGrid
              DoseInGrid<-!is.na(match(dose,data@doseGrid))
              if (DoseInGrid==TRUE){
                ##Find which dose is this in the dose Grid
                EIx<-which(dose==data@doseGrid)
                ##Return corresponding expected efficacy values from mcmc samples
                return(Effsamples@data$ExpEff[,EIx])} else {##if dose not in doseGrid do linear Interploation
                  ## check if this dose is within doseGrid
                  stopifnot(dose <= max(data@doseGrid), dose >= min(data@doseGrid))
                  Ixd0 <- max(which((dose-data@doseGrid) > 0))
                  Ixd1<-min(which((dose-data@doseGrid) < 0))
                  ExpEffd0<-Effsamples@data$ExpEff[,Ixd0]
                  ExpEffd1<-Effsamples@data$ExpEff[,Ixd1]
                  InterpolExpEff<- ExpEffd0+(ExpEffd1-ExpEffd0)*((dose-data@doseGrid[Ixd0])/(data@doseGrid[Ixd1]-data@doseGrid[Ixd0]))
                  return(InterpolExpEff)}
            },
            useFixed=useFixed,
            useRW1=useRW1,
            designW=designW,
            RWmat=RWmat,
            RWmatRank=RWmatRank)}
## ---------------------------------------------------------------------------------------------------------

