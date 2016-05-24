#' Draw values from multinomial distribution
#' 
#' Draws values from a mutinomial distribution and matches them with a row in the pmatrix (attribute profile matrix) 
#' to obtain the alpha values
#' @param pmatrix a matrix indicating all possible attribute profiles
#' @param y number of random vectors to draw
#' @param nus0 a numeric vector of length \code{nclasses} for nu parameter estimates
#' @return a numeric vector of alpha values
DrawAlphas <- function(pmatrix, y, nus0){
  pmatrix[which(rmultinom(y, 1, nus0) ==1),]
}
#' Draw latent variable values
#' 
#' Produces values of alphas(s) (i.e., latent variable(s)) by using a rejection/acceptance Metropolis-Hastings procedure.
#' @param nattributes a numeric value indicating number of alphas to be drawn
#' @param class0 a numeric vector of previous classes
#' @param nus0 a numeric vector of length \code{nclasses} for nu parameter estimates
#' @param observations a data frame or matrix of dichotomous responses
#' @param nobservations a numeric value indicating number of rows of the observation data frame or matrix
#' @param threshold.values an nclasses by nitems numeric matrix with appropriate item threshold values
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @return a numeric vector of classes indicating one of the possible attribute profile for each student
DrawClasses <- function(nattributes, class0, nus0, observations, nobservations, threshold.values, pmatrix) {
  alphas <- t(sapply(rep(1, nobservations), function(y) DrawAlphas(pmatrix, y, nus0)))
  class1 <- apply(alphas, 1, function(y) which(apply(pmatrix, 1, function(x) all(x == y))))
  probs0 <- InvLogit(threshold.values[class0,])
  probs0 <- ifelse(observations == 1, probs0, 1 - probs0)
  dcm0 <- apply(probs0, 1, function(x) prod(x, na.rm=T))
  probs1 <- InvLogit(threshold.values[class1,])
  probs1 <- ifelse(observations == 1, probs1, 1 - probs1)
  dcm1 <- apply(probs1, 1, function(x) prod(x, na.rm=T))
  accept <- ifelse(dcm0==0,1,dcm1/dcm0)
  accept <- ifelse(accept>1,1,accept)
  accept <- ifelse(runif(nobservations)<accept,1,0) 
  class1 <- ifelse(accept,class1,class0) 
  return(class1)
}

#' Perform one iteration of MCMC procedure
#'
#' If applicable, randomly samples new set of parameter estimates, obtains applicable estimates
#' and uses those to calculate threshold values for both items and latent variables, draws
#' new set of alpha values.
#' @param nattributes numberic value for number of attributes
#' @param class0 The previous value of attribute profile for each respondent
#' @param estimates0 a numeric vector of parameter estimates
#' @param threshold.labels an nclasses by nitems character matrix with appropriate threshold labels
#' @param lambda.equations equations for lambda parameters
#' @param is.pi.r If \code{FALSE} (the default), parameter values are the type of \code{taus} and \code{nus} or 
#' \code{lambdas} and \code{gammas} else they are the type \code{pis} and \code{rs} as used in NC-RUM parameterization
#' @param parameter.means a numerical vector of calibrated item and structural parameters
#' @param parameter.acov a numerical matrix of covariances of item and structural parameters
#' @param observations a data frame or matrix of dichotomous responses
#' @param nobservations a numeric value of number of observations
#' @param is.parameter.randomized if true parameter estimates are randomized using acov matrix
#' @param qmatrix a data frame or matrix of 1s and 0s indicating relation between items and attributes. 
#' This matrix specifies which items are required for mastery of each attribute (i.e., latent variable). 
#' A matrix must be a size of \code{nItems X nAttributes}
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @return a list of newly sampled classes and parameter estimates
iterate <-function(nattributes, class0, estimates0, threshold.labels, lambda.equations
                   , is.pi.r, parameter.means, parameter.acov, observations, nobservations
                   , is.parameter.randomized, qmatrix, pmatrix) {

  if (is.parameter.randomized)
    estimates0 <- SampleParameterEstimates(parameter.means, parameter.acov)
  if (!is.null(lambda.equations)){
    threshold.info <- GetThresholdValuesKernel(lambda.equations, estimates0, threshold.labels)
    taus <- -threshold.info$taus
    threshold.values <- -threshold.info$threshold.values
    mus <- GetMusFromGammas(estimates0, pmatrix)
    estimates0 <- c(taus, mus)
  }
  if (is.pi.r){
    taus <- GetThresholdValuesKernelPiR(estimates0, qmatrix, pmatrix, threshold.labels)
    mus <- GetMusFromGammas(estimates0, pmatrix)
    estimates0 <- c(taus, mus)
  }
  parms <- InitializeParameters(estimates0, nrow(pmatrix))
  taus0 <- parms[['taus']]
  nus0 <- parms[['nus']]
  if (is.null(lambda.equations)){
    threshold.values <- GetThresholdValues(threshold.labels, taus0)
  }
  class1 <- DrawClasses(nattributes, class0, nus0, observations, nobservations, threshold.values, pmatrix)
  return(list(new.class = class1, estimates = estimates0))
}

#' Randomly sample parameter estimates
#' 
#' Given parameter estimates and asymptotic covaraince matrix of parameter estimates, randomly samples a new
#' set of estimates
#' @param means a numerical vector of calibrated item and structural parameters
#' @param acov a numerical matrix of covariances of item and structural parameters
#' @return A vector of randomly sampled parameter estimate values
SampleParameterEstimates <- function(means, acov){
  #TODO: check for non-singularity of acov matrix
  if (det(acov) != 0)
    L= chol(acov)
  else
    L = t(ldl(acov))    #chol decomposition of acov (L*L` = acov solve for L, square root of acov matrix)
  Z = rnorm(matrix(0, 1, length(means)))		#random standard normal variates (theta ~ N(means, acov) 
                                              # where theta is vector of all parameters- taus and mus)
  Z = matrix(Z, 1, length(means))			#placing in matrix form
  simulated.values = c(Z%*%L) +  means#getting simulated values (theta = N(0,1)*sqrt(acov) + means)
  return(simulated.values)
}

#' Initialize parameter estimates
#' 
#' Obtains appropriate parameter estimates for a non-kernel model (e.g., taus and nus)
#' @param estimates a vector of parameter estimates
#' @param nclasses a numeric value representing the number of unique attribute profiles 
InitializeParameters <- function(estimates, nclasses){
  taus <- unlist(estimates[grep('^tau', names(estimates))])
  mus <- estimates[grep('^mu', names(estimates))]
  if (length(mus) < nclasses)
    mus <- c(mus, rep(0, nclasses -length(mus)))
  nus <- GetClassProbsFromMus(mus)

  list("taus" = taus, "nus" = nus)
}

#' Performs MCMC routine for DCM
#' 
#' @param observations a data frame or matrix of dichotomous responses
#' @param nattributes numeric value of number of attributes
#' @param qmatrix a data frame or matrix of 1s and 0s indicating relation between items and attributes. 
#' This matrix specifies which items are required for mastery of each attribute (i.e., latent variable). 
#' A matrix must be a size of \code{nItems X nAttributes}
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @param parameter.means a numerical vector of calibrated item and structural parameters
#' @param parameter.acov a numerical matrix of covariances of item and structural parameters
#' @param nobservations a numeric value indicating number of rows of the observation data frame or matrix
#' @param nreps The number of iterations in MCMC per chain
#' @param initial.class The initial value of attribute profile for each respondent
#' @param nchains The number of chains in MCMC
#' @param threshold.labels an nclasses by nitems character matrix with appropriate item threshold labels
#' @param lambda.equations lambda parameter equations
#' @param is.pi.r If \code{FALSE} (the default), parameter values are the type of \code{taus} and \code{nus} or 
#' \code{lambdas} and \code{gammas} else they are the type \code{pis} and \code{rs} as used in NC-RUM parameterization
#' @param is.parameter.randomized if true parameter estimates are randomized using acov matrix
#' @param parameterization.method optional character string of parameterization method used to calibrate parameters
#' @param percent.reps.to.discard The percent of iterations to be discarded
#' @return a list of class and parameter data frame containing all accepted iteraction of MCMC
mcmc <- function(observations, nattributes, qmatrix, pmatrix, parameter.means
                 , parameter.acov, nobservations, nreps, initial.class, nchains
                 , threshold.labels, lambda.equations, is.pi.r
                 , is.parameter.randomized, parameterization.method, percent.reps.to.discard){
                   

  chain.class.results <- list()
  chain.parameter.results <- list()
  cols.to.discard = ceiling(nreps*percent.reps.to.discard*0.01)
  # begin chains
  for (chain in 1:nchains){
    # sample from the parameter estimates (mean vector) and acov matrix (cov matrix)
    estimates <- parameter.means
    class.out <- matrix(NA, nobservations, nreps)
    class.out[, 1] <- initial.class
    class0 <- initial.class
    estimates.out <- matrix(NA, length(estimates), nreps)
    rownames(estimates.out) <- names(estimates)
    estimates.out[, 1] <- estimates
    
    # begin replications 
    for(i in 2:nreps){
      # perform one iteration
      out <- iterate(nattributes, class0, estimates, threshold.labels, lambda.equations, is.pi.r
                     , parameter.means, parameter.acov, observations, nobservations
                     ,is.parameter.randomized, qmatrix, pmatrix)
      class.out[, i] <- class0 <- out[["new.class"]]
      estimates.out[, i] <- out[["estimates"]]
    }
    chain.class.results[[chain]] <- class.out[, -(1:cols.to.discard)]
    chain.parameter.results[[chain]] <- estimates.out[, -(1:cols.to.discard)]
  }
  names(chain.class.results) <- names(chain.parameter.results) <- 1:nchains
  return(list(class.result = chain.class.results, parameter.result = chain.parameter.results))
}

#' Calculate attribute probabilities
#' 
#' Calculates the attribute (i.e., latent variable) probabilities for all students
#' @param class.result a list of classes containing all interactions and chains of MCMC
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
GetAttributesProbabilitiesMCMC <- function(class.result, pmatrix){
  #theta_values <- ddply(narrow.results, .(rep, id, group), print(group))
  narrow.results <- with(class.result, LongFormatResults(class.result, 'klass'))
  alpha.values <- with(narrow.results, ddply(narrow.results, .(chain, rep, id, klass), function(x) pmatrix[x$klass,]))
  narrow.alpha.values <- with(alpha.values, melt(alpha.values, id.vars = c("chain", "rep", "id", "klass"), variable.name = "attr.name", value.name = "attr.value"))
  narrow.alpha.results <- with(narrow.alpha.values, ddply(narrow.alpha.values
                                , .(id, attr.name), summarize, means = mean(attr.value), sds = sd(attr.value)))
  narrow.alpha.results.melted <- with(narrow.alpha.results, melt(narrow.alpha.results
                                      , id.vars = c("id", "attr.name"), variable.name = "moment", value.name = "moment.value"))
  alpha.results <- with(narrow.alpha.results.melted, dcast(narrow.alpha.results.melted, id ~ attr.name + moment, value.var = 'moment.value'))
  return(alpha.results)
}

#' Calculate attribute profile probabilities
#' 
#' Calculates the attribute profile probabilities for all students
#' @param class.result a list of classes containing all interactions and chains of MCMC
#' @param nclasses a numeric value representing the number of unique attribute profiles 
GetClassProbabilitiesMCMC <- function(class.result, nclasses){
  narrow.results <- with(class.result, LongFormatResults(class.result, 'klass'))
  narrow.class.probs <- with(narrow.results, ddply(narrow.results, .(chain, id, klass), summarize, prob = length(klass)/max(rep)))
  narrow.class.probs <- with(narrow.class.probs, ddply(narrow.class.probs, .(id, klass), summarize, prob = mean(prob)))
  narrow.class.probs <- with(narrow.class.probs, rbind(narrow.class.probs, data.frame(id = rep(0, nclasses), klass = 1:nclasses, prob = rep(0, nclasses))))
  max_ <- with(narrow.class.probs, ddply(narrow.class.probs, .(id), summarize, max.class = klass[which(max(prob, na.rm=T) == prob)]))
  class.probs <- with(narrow.class.probs, dcast(narrow.class.probs, id ~ klass, value.var = 'prob', fill = 0))
  class.probs.max <- merge(class.probs, max_, by = "id")
  return(class.probs.max[!class.probs.max$id == 0,])
}

#' Calculate results for attributes and attribute profiles
#' 
#' Calculates results for attributes and attribute profiles for all students based on MCMC output
#' @param parameter.result a list parameter estimates for each iterations and chains
GetParameterResultsMCMC <- function(parameter.result){
  parameter.result <- with(parameter.result, LongFormatResults(parameter.result, 'estimate'))
  parameter.result <- with(parameter.result, ddply(parameter.result, .(id), summarize, parm.means = mean(estimate), parm.sds = sd(estimate)))
  return(parameter.result)
}

#' Tranform dataframe to long format
#' 
#' Given a dataframe in wide format, transform to long format
#' @param d a list of MCMC output either classes or parameter estimates
#' @param result.type optional character strings containing result.type: \code{class} or \code{parameter}
LongFormatResults <- function(d, result.type = "class"){
  long.d <- ldply(d, function(x) melt(x, value.name = result.type))
  long.d <- rename(long.d, c(".id" = "chain", "Var1" = "id", "Var2" = "rep"))
  return(long.d)
}

#' Score dichotomous response data using DCM and MCMC
#'    
#' Given dichotomous response data, a Q-matrix (relation between items and attributes) and calibrated item and structural parameters
#' \code{ScoreDCM} estimates posterior probabilities of attribute profiles of respondents using a Diagnostic Classification Model (DCM)
#' and Markov Chain Monte Carlo (MCMC) method. The estimation procedure uses the loglinear cognitive diagnostic modeling (LDCM)
#' framework that allows for the estimation of a host of DCMs such as DCM, DINA, C-RUM, NIDO, NIDA, NC-RUM etc.
#' @param observations a data frame or matrix of dichotomous responses in wide format (\code{nobservations X nitems})
#' @param qmatrix a data frame or matrix of 1s and 0s indicating relation between items and attributes. This matrix specifies which items are required
#'  for mastery of each attribute (i.e., latent variable). The matrix must be a size of \code{nitems X nattributes}
#' @param parameter.means a numerical vector of calibrated item and structural parameters. Values must be in the order of 
#' \code{\link{GetParameterNames}} if parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values 
#' are used must be in order of \code{\link{GetKernelParameterNames}}
#' @param parameter.acov optional matrix of covariances of all model parameters. If \code{NULL} (the default) model parameters are not randomized 
#' for each iteration of MCMC
#' @param parameterization.method optional character string of parameterization method used to calibrate parameters. 
#' If not specified then the default will be set to \code{Mplus}
#' @param is.kernel.parameters If \code{FALSE} (the default), parameter values are of type \code{taus} and \code{nus} else they are of type 
#' kernel parameters, i.e., \code{lambdas} and \code{gammas}
#' @param model.type If \code{is.kernel.parameter} is \code{TRUE}, model type must be one of \code{DCM}, \code{DINA}, \code{CRUM}, \code{DINO},
#' \code{NIDO}, \code{NCRUM}. Kernel parameters are different for each model type
#' @param nreps The number of iterations in MCMC per chain. If not specified the default is 1000
#' @param nchains The number of chains in MCMC. If not specified the default is 1
#' @param initial.class The initial value of the attribute profile for each respondent. If not specified the default is 1
#' @param percent.reps.to.discard The percent of iterations to be discarded. If not specifed then the default is 5
#' 
#' @return ScoreDCM returns an object of class \code{\link{dcm.scorer.class}}; a list consisting of
#' \describe{
#'  \item{inputs}{a list of all input arguments to the funciton like observations, Q-matrix, parameter.means, etc.}
#'  \item{mcmc.inputs}{a list of all mcmc inputs}
#'  \item{results}{an object of \code{\link{all.results.class}} class; a list consisting of
#'    \describe{
#'      \item{attribute.profile.result}{an object of \code{\link{attribute.profile.class}}}
#'      \item{attribute.result}{an object of \code{\link{attribute.class}}}
#'      \item{parameter.result}{an object of \code{\link{parameter.class}}}
#'    }
#'  }
#'  \item{mcmc.outputs}{a list of all mcmc output}
#' }
#' The function \code{\link{summary}} is used to obtain and print a summary of
#' MCMC runs in the form of probabilities of mastering each attribute and attribute profile probabilities.
#' The function \code{\link{plot}} is used to plot the aggregated mean of both attribute mastery (type = "attr.means")
#' and attribute profile probabily (type = "attr.profile.means") across all respondents. Other plot options include
#' attribute mastery (type = "attr.profiles") and attribute profile probabilities (type = "attr.profile.profiles") 
#' of individual respondents. 
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
#' @usage ScoreDCM(observations, qmatrix, parameter.means, parameter.acov= NULL, 
#' parameterization.method = "Mplus", is.kernel.parameters = FALSE, model.type = NULL,
#' nreps = 1000 , nchains = 1, initial.class = 1,  percent.reps.to.discard = 5)
#' @examples
#' 
#' \dontrun{
#' 
#' data(obervations.test)
#' data(qmatrix.test)
#' parameter.names <- GetParameterNames(qmatrix.test, ncol(qmatrix.test))
#' parameter.names
#' data(parameter.means.DCM.Mplus.test)
#' obj <- ScoreDCM(observations = observations.test, qmatrix = qmatrix.test
#' , parameter.means = parameter.means.DCM.Mplus.test)
#' summary(obj)
#' plot(obj)
#' }
#' @export
ScoreDCM <- function(observations, qmatrix, parameter.means, parameter.acov= NULL, 
                     parameterization.method = "Mplus", is.kernel.parameters = FALSE, model.type = NULL,
                     nreps = 1000 , nchains = 1, initial.class = 1,  percent.reps.to.discard = 5){

  kAllModelTypes <- c('DCM', 'DINA', 'CRUM', 'DINO', 'NIDO', 'NCRUM')
  is.pi.r <- FALSE
  if (!hasArg(observations))
    stop("argument \"observations\" is missing, with no default")
  if (!hasArg(qmatrix))
    stop("argument \"qmatrix\" is missing, with no default")
  if (!hasArg(parameter.means))
    stop("argument \"parameter.means\" is missing, with no default")
  if (is.null(observations) | (!(is.data.frame(observations) | is.matrix(observations))))
    stop("argument \"observations\" must be data frame or a matrix")
  if (nrow(observations) == 0)
    stop("argument \"observations\" can not be empty")
  if (!all(sapply(observations, is.numeric)))
    stop("argument \"observations\" must be a matrix of all numeric values")
  if (all(sapply(observations, function(x) x!=1 & x!=0)))
    stop("argument \"observations\" must be a matrix of 0s and 1s")
  if (!is.matrix(qmatrix))
    stop("argument \"qmatrix\" must be a matrix")
  if (nrow(qmatrix) == 0 | ncol(qmatrix) == 0)
    stop("argument \"qmatrix\" should not be empty")
  if (!all(sapply(qmatrix, is.numeric)))
    stop("argument \"qmatrix\" must be a matrix of all numeric values")
  if (all(sapply(qmatrix, function(x) x!=1 & x!=0)))
    stop("argument \"qmatrix\" must be a matrix of 0s and 1s")
  if (nrow(qmatrix) != ncol(observations))
    stop(sprintf("argument \"qmatrix\" must be a matrix %d rows", ncol(observations)))
  if (!all(is.numeric(parameter.means)))
    stop("parameter.means must be vector of numeric values")
  if (!all(sapply(qmatrix, is.numeric)))
    stop("parameter.acov must be matrix of numeric values")
  if (is.kernel.parameters){
    if (is.null(model.type))
      stop(sprintf("model.type must be one of: %s", paste(kAllModelTypes, collapse = ', ')))
    if(!(model.type %in% kAllModelTypes))
      stop(sprintf("model.type must be one of: %s", paste(kAllModelTypes, collapse = ', ')))
  }
  nobservations <- nrow(observations)
  nitems <- nrow(qmatrix)
  nattributes <- ncol(qmatrix)
  nclasses <- 2^nattributes
  
  if (!is.null(rownames(qmatrix))){
    item.names <- rownames(qmatrix)
  } else{
    item.names <- sprintf('item%d', 1:nitems)
  }
  if (!is.null(colnames(qmatrix))){
    attribute.names <- colnames(qmatrix)
  } else{
    attribute.names <- sprintf('attribute%d', 1:nattributes)
  }
  dimnames(qmatrix) <- list(item.names, attribute.names)
  pmatrix <- GetAttributeProfiles(nattributes, attribute.names)
  #required parameters for a given model type
  if (!is.kernel.parameters){
    required.parms <- GetParameterNames(qmatrix, nattributes, parameterization.method)
    names(parameter.means) <- required.parms
    if (parameterization.method == "Mplus"){
      parameter.means[grep('^tau', names(parameter.means))] <- -parameter.means[grep('^tau', names(parameter.means))]
    }
    lambda.equations <- NULL
  } else {
    kernel.parms <- GetKernelParameterNames(qmatrix, nattributes, model.type)
    if (model.type == "NCRUM"){
      names(parameter.means) <- kernel.parms$pis.rias
      is.pi.r <- TRUE
    } else {
      names(parameter.means) <- c(kernel.parms$lambdas$lambda.name, kernel.parms$gammas)
      lambda.equations <- kernel.parms$lambdas$lambda.equations
    }
  }
  if(length(parameter.means) != length(required.parms))
    stop(sprintf("parameter means should be a length of %d", length(required.parms)))
  if(!is.null(parameter.acov)){
    is.parameter.randomized = TRUE
    if ((nrow(parameter.acov) != ncol(parameter.acov)) | (ncol(parameter.acov) != length(required.parms)))
      stop(sprintf("parameter acov should be a matrix of [%d, %d]", length(required.parms), length(required.parms)))
  } else{
    is.parameter.randomized = FALSE
  }
  initial.class <- rep(initial.class, nobservations)
  threshold.labels <- GetThresholdLabels(qmatrix, pmatrix)
  mcmc.out <- mcmc(observations, nattributes, qmatrix, pmatrix, parameter.means
                   , parameter.acov, nobservations, nreps, initial.class, nchains
                   , threshold.labels, lambda.equations, is.pi.r, is.parameter.randomized
                   , parameterization.method, percent.reps.to.discard)
  attribute.profile.out <- attribute.profile.class(results = GetClassProbabilitiesMCMC(mcmc.out$class.result, nclasses), attribute.profile.matrix = pmatrix)
  attribute.out <- attribute.class(results = GetAttributesProbabilitiesMCMC(mcmc.out$class.result, pmatrix))
  parameter.out <- parameter.class(results = GetParameterResultsMCMC(mcmc.out$parameter.result))
  inputs <- list(observations = observations, qmatrix = qmatrix, parameter.means = parameter.means,
                 nobservations = nobservations, nitems = nitems, nattributes = nattributes, nclasses = nclasses,
                 attribute.profiles = pmatrix, parameterization.method = parameterization.method,
                 is.parameter.randomized = is.parameter.randomized, parameter.acov = parameter.acov)
  mcmc.inputs <- list(nreps = nreps, nchains = nchains, burn = percent.reps.to.discard,
                      initial.class = initial.class)
  outputs <- all.results.class(attribute.profile.result = attribute.profile.out , attribute.result = attribute.out,
                  parameter.result= parameter.out)
  score.dcm <- dcm.scorer.class(inputs = inputs, mcmc.inputs = mcmc.inputs, results = outputs, mcmc.output = mcmc.out)                
  return(score.dcm)
}

