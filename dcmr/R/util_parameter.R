#' Attribute Profiles
#'
#' Given a number of attributes, it generates all possible attribute profiles (response patterns of latent variables)
#' @param nattributes a numeric value for number of attributes
#' @param attribute.names optional character vector of attribute names
#' @param profile.names optional character vector of attribute profile names
#' @return a matrix of dimension \code{2^nattributes X nattributes} containing binary values
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
#' @usage GetAttributeProfiles(nattributes, attribute.names = NULL, profile.names = NULL)
#' @export 
#' @examples
#'
#' \dontrun{
#'
#'  GetAttributeProfiles(nattributes=3)
#' }
#'
GetAttributeProfiles <- function(nattributes, attribute.names = NULL, profile.names = NULL){
  nclasses <- 2^nattributes
  numeric.to.binary <- as.numeric(unlist(strsplit(as.character(binary(0:(nclasses-1))), "")))
  pmatrix <- matrix(numeric.to.binary, nclasses, nattributes, byrow = T, dimnames = list(profile.names, attribute.names))
  pmatrix
}

#' Parameter Names for All DCM Models
#' 
#' Given a Q-matrix and parameterization method it generates item and structural parameter names. These are non-kernel parameters: 
#' item thresholds (\code{taus}) and latent variable thresholds (\code{nus})
#' @param qmatrix a matrix of size (nitems X nattributes) that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param nattributes a numeric value of number of attriutes
#' @param parameterization.method optional character string of parameterization method used to calibrated parameters. 
#' If not specified then the default will be set to  \code{Mplus}
#' @return a string vector of threshold parameter names
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
#' @usage GetParameterNames(qmatrix, nattributes, parameterization.method = 'Mplus')
#' @export
#' @examples
#' 
#' \dontrun{
#' 
#'  GetParameterNames(qmatrix = qmatrix.test, nattributes=3)
#' }
GetParameterNames <- function(qmatrix, nattributes, parameterization.method = 'Mplus'){
  pmatrix <- GetAttributeProfiles(nattributes)
  threshold.out <- NULL
  for (c in 1:nrow(pmatrix)){
    for (i in 1:nrow(qmatrix)){
      # numeric vector of required attributes for an item i given qmatrix and class c 
      required.attributes <-  which(qmatrix[i,] == 1)[which(qmatrix[i,] == 1) %in% which(pmatrix[c,] == 1)]    
      if (length(required.attributes) == 0) required.attributes = "0"
      kernel <- paste(required.attributes, collapse = "")
      threshold.out <- c(threshold.out, sprintf("tau%d_%s", i, kernel))
    }
  }
  if(parameterization.method == "Mplus")
    mus.out <- sprintf("mu_%d", 1:(nrow(pmatrix) - 1))
  else 
    mus.out <- sprintf("mu_%d", 1:(nrow(pmatrix)))
  return(c(unique(threshold.out), unique(mus.out)))
}

#' Obtain name of a particular lambda parameter
#' 
#' Gets the name of a particular lambda parameter and produces a string of the
#' form: l_item_effect_attribtues.
#' @param item a numeric value representing item number
#' @param attributes a String of the attributes associated with this gamma
#' @return a String of the form: l_item_effect_attribtues.
GetLambdaName <- function(item, attributes){
  return(sprintf("l%d_%d_%s", item, nchar(attributes), paste(attributes, collapse = "")))
}

#' Gets lambda names for a given item
#' 
#' @param item.index.intercept a numeric value representing item number for the intercept lambda
#' @param item.index.main.effect a numeric value representing item number for the main effects lambda
#' @param required.attributes a numeric vector of attributes required to have mastered to get item correct
#' @return a String vector of lambda names for one item
GetLambdaNamesForItem <- function(item.index.intercept, item.index.main.effect, required.attributes){
  final.lambdas <- sprintf("l%d_0", item.index.intercept)
  if (any(grepl("&", required.attributes))){
    required.attributes <- as.numeric(regmatches(required.attributes, regexpr('[0-9]+', required.attributes)))
    for (i in 1:length(required.attributes)){
      attrs <- combn(required.attributes, i, function(x) paste(x, collapse=""))
      tmp.lambda <- sapply(attrs, function(x) GetLambdaName(item.index.main.effect, x))
      final.lambdas <- c(final.lambdas, tmp.lambda)
    }
  } else if(any(grepl("\\|", required.attributes))){
    required.attributes <- as.numeric(regmatches(required.attributes, regexpr('[0-9]+', required.attributes)))
    tmp.lambda <- sapply(required.attributes, function(x) GetLambdaName(item.index.main.effect, x))
    final.lambdas <- c(final.lambdas, tmp.lambda)
  } else if(length(required.attributes) == 1) {
    final.lambdas <- c(final.lambdas, sprintf("l%d_1_%s", item.index.main.effect,  required.attributes))
  }
  return(final.lambdas)
}
#' Gets names of lambda parameters
#' 
#' @param item.indexes.intercepts a numeric vector indicating which item is assoicated with the lamdba intercept
#' at the given index and thus allowing for equality constraints. For example, with nitems = 3,
#' this is c(1, 1, 1), that implies that intercept lambdas across items should be contrained to equality.
#' @param item.indexes.main.effects a numeric vector indicating which item is assoicated with the lamdba main effect
#' at the given index and thus allowing for equality constraints. For example, with nitems = 3,
#' this is c(1, 1, 1), that implies that main effect (and higher order) lambdas 
#' across items should be contrained to equality.
#' @param required.attributes.lambda a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to lambda parameters.
#' @return a String vector of names of lambda parameters
GetLambdaNames <- function(item.indexes.intercepts, item.indexes.main.effects, required.attributes.lambda){
  nitems <- ncol(required.attributes.lambda)
  nclasses <- nrow(required.attributes.lambda)
  all.lambdas <- NULL
  lambda.equations <- matrix(NA, nclasses, nitems) 
  for (i in 1:ncol(required.attributes.lambda)){
    lambdas <- NULL
    for (c in 1:nrow(required.attributes.lambda)){
      # numeric vector of required attributes for an item i given qmatrix and class c 
      required.attributes <- unlist(strsplit(required.attributes.lambda[c, i], ""))
      cur.lambdas <- GetLambdaNamesForItem(item.indexes.intercepts[i]
                                           , item.indexes.main.effects[i]
                                           , required.attributes)
      lambdas <- c(lambdas, cur.lambdas)
      lambda.equations[c, i] <- paste(cur.lambdas, collapse = " + ")
    }
    all.lambdas <- c(all.lambdas, lambdas)
  }
  return(list(lambda.names = unique(all.lambdas), lambda.equations = lambda.equations))
}

#' Generate required attributes
#' 
#' Generates a matrix of required attributes
#' @param qmatrix a nitems by nattributes matrix that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @return a nclasses by nitems character matrix with Strings representing the required attributes for a given item and class.
GetRequiredAttributes <- function(qmatrix, pmatrix){
  apply(qmatrix, 1, function(y) {
    aaply(pmatrix, 1, function(x) 
      paste(which(y ==1)[which(y== 1) %in% which(x== 1)], collapse = ""))}
  )
}

#' Calculate required attributes for lamdbas fully specified DCM
#' 
#' Calculates required attributes for a model parameterized in the LCDM framework using the fully
#' specified DCM model.
#' @param required.attributes.tau a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to tau parameters.
#' @return a nclasses by nitems matrix with each cell containing a String representing the
#' attributes required for the given class and item specific to lambda parameters for the a fully
#' specific DCM model.
GetRequiredAttributesLambdaFullDCM <- function(required.attributes.tau){
  apply(required.attributes.tau, c(1,2), function(x) ifelse(nchar(x) > 1, paste(unlist(strsplit(x, "")), collapse="&"), x))
}
#' Calculate required attributes for lambdas for a DINA model.
#' 
#' Calculates required attributes for a model parameterized in the LCDM framework using the DINA
#' model
#' @param required.attributes.tau a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to tau parameters.
#' @return a nclasses by nitems matrix with each cell containing a String representing the
#' attributes required for the given class and item specific to lambda parameters for the DINA model.
GetRequiredAttributesLambdaDina <- function(required.attributes.tau){
  apply(required.attributes.tau, 2, function(x) ifelse(nchar(x) == max(nchar(x)), "e", ""))
}
#' Calculate required attributes for lambdas for a DINO model.
#' 
#' Calculates required attributes for a model parameterized in the LCDM framework using the DINO
#' model
#' @param required.attributes.tau a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to tau parameters.
#' @return a nclasses by nitems matrix with each cell containing a String representing the
#' attributes required for the given class and item specific to lambda parameters for the DINO model.
GetRequiredAttributesLambdaDino <- function(required.attributes.tau){
  apply(required.attributes.tau, 2, function(x) ifelse(x == "", "", "e"))
}
#' Calculate required attributes for lambdas for a CRUM model.
#' 
#' Calculates required attributes for a model parameterized in the LCDM framework using the CRUM
#' model
#' @param required.attributes.tau a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to tau parameters.
#' @return a nclasses by nitems matrix with each cell containing a String representing the
#' attributes required for the given class and item specific to lambda parameters for the CRUM model.
GetRequiredAttributesLambdaCrum <- function(required.attributes.tau){
  apply(required.attributes.tau, c(1,2), function(x) ifelse(nchar(x) > 1, paste(unlist(strsplit(x, "")), collapse="|"), x))
}
#' Calculate required attributes for lambdas for a NIDO model.
#' 
#' Calculates required attributes for a model parameterized in the LCDM framework using the NIDO
#' model
#' @param required.attributes.tau a nclasses by nitems matrix with each cell containing a String
#' representing the attributes required for the given class and item specific to tau parameters.
#' @return a nclasses by nitems matrix with each cell containing a String representing the
#' attributes required for the given class and item specific to lambda parameters for the NIDO model.
GetRequiredAttributesLambdaNido <- function(required.attributes.tau){
  apply(required.attributes.tau, c(1,2), function(x) ifelse(nchar(x) > 1, paste(unlist(strsplit(x, "")), collapse="|"), x))
}

#' Obtain name of a particular gamma parameter
#' 
#' Produces a string of the form: g_effect_attribtues representing gamma parameters
#' @param attributes a string of the attributes associated with this gamma parameter 
#' @return a string of the form g_effect_attribtues.
GetGammaName <- function(attributes){
  return(sprintf("g_%d_%s", nchar(attributes), paste(strsplit(attributes, "")[[1]], collapse = "")))
}

#' Obtain names of all gamma parameters
#' 
#' Gets names of all gamma parameters
#' @param required.attributes a numeric vector of attributes required to have mastered to get item correct
#' @return a String vector of names for gamma parameters
GetGammaNames <- function(required.attributes){
  final.gammas <- NULL
  if (length(required.attributes) > 1){
    for (i in 1:length(required.attributes)){
      attrs <- combn(required.attributes, i, function(x) paste(x, collapse=""))
      tmp.gamma <- sapply(attrs, function(x) GetGammaName(x)) 
      final.gammas <- c(final.gammas, tmp.gamma)
    }
  }
  if (length(required.attributes) == 1) {
    final.gammas <- c(final.gammas, sprintf("g_1_%d", required.attributes))
  }
  return(unique(final.gammas))
}

#' Kernel Parameter Names for all DCM Models
#' Given a Q-matrix and model type it generates item and structural parameter names. These are kernel parameters: 
#' item thresholds (\code{lambdas}) and latent variable thresholds (\code{gammas})
#' @param qmatrix a matrix of size (nitems X nattributes) that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param nattributes a numeric value of number of attriutes
#' @param model.type a string value of \code{DCM}, \code{DINA}, \code{CRUM}, \code{DINO},
#' \code{NIDO}, \code{NCRUM}.
#' If not specified then the default will be set to  \code{DCM}
#' @return a string vector of threshold parameter names
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
#' @usage GetKernelParameterNames(qmatrix, nattributes, model.type = "DCM")
#' @export
#' @examples
#' 
#' \dontrun{
#' 
#' GetKernelParameterNames(qmatrix = qmatrix.test, nattributes=3)
#' }
#' 
GetKernelParameterNames <- function(qmatrix, nattributes, model.type = "DCM"){
  pmatrix <- GetAttributeProfiles(nattributes)
  nitems <- nrow(qmatrix)
  nclasses <- nrow(pmatrix)
  gamma.names <- GetGammaNames(which(pmatrix[nclasses,] == 1))
  required.attributes.tau <- GetRequiredAttributes(qmatrix, pmatrix)
  if (model.type == "DCM"){
    required.attributes.lambda <- GetRequiredAttributesLambdaFullDCM(required.attributes.tau)
    item.indexes.intercepts <- c(1:nitems)
    item.indexes.main.effects <- c(1:nitems)
  }
  if (model.type == "DINA"){
    required.attributes.lambda <- GetRequiredAttributesLambdaDina(required.attributes.tau)
    item.indexes.intercepts <- c(1:nitems)
    item.indexes.main.effects <- c(1:nitems)
  }
  if (model.type == "DINO"){
    required.attributes.lambda <- GetRequiredAttributesLambdaDino(required.attributes.tau)
    item.indexes.intercepts <- c(1:nitems)
    item.indexes.main.effects <- c(1:nitems)
  }
  if (model.type == "CRUM"){
    required.attributes.lambda <- GetRequiredAttributesLambdaCrum(required.attributes.tau)
    item.indexes.intercepts <- c(1:nitems)
    item.indexes.main.effects <- c(1:nitems)
  }
  if (model.type == "NIDO"){
    required.attributes.lambda <- GetRequiredAttributesLambdaNido(required.attributes.tau)
    item.indexes.intercepts <- c(1:nitems)
    item.indexes.main.effects <- rep(1, nitems)
  }
  if (model.type == "NCRUM"){
    pi.i.stars <- sprintf("pi.%d.star", 1:nrow(qmatrix))
    r.ia.stars <- unlist(sapply(1:nrow(qmatrix), function(x) sprintf("r.%d.%d.star", x, which(qmatrix[x, ] == 1))))
    return(list(pis.rias = c(pi.i.stars, r.ia.stars), gammas = gamma.names))
  } else {
    lambda.names.structure <- GetLambdaNames(item.indexes.intercepts
                                             , item.indexes.main.effects
                                             , required.attributes.lambda)
  }
  return(list(lambdas = lambda.names.structure, gammas = gamma.names))
}

#' Obtain item threshold labels
#' 
#' Gets item threshold labels for a given qmatrix and pmatrix
#' @param qmatrix a numeric nitems by nattributes matrix that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @return an nclasses by nitems character matrix with appropriate threshold labels.
#' @export
GetThresholdLabels <- function(qmatrix, pmatrix){
  threshold.labels <- NULL
  for (i in 1:nrow(qmatrix)){
    item.threshold.labels <- NULL
    for (c in 1:nrow(pmatrix)){
      # numeric vector of required attributes for an item i given qmatrix and class c 
      required.attributes <-  which(qmatrix[i,] == 1)[which(qmatrix[i,] == 1) %in% which(pmatrix[c,] == 1)]    
      if (length(required.attributes) == 0) required.attributes = "0"
      item.threshold.labels <- c(item.threshold.labels, sprintf("tau%d_%s", i, paste(required.attributes, collapse = "")))
    }
    threshold.labels <- cbind(threshold.labels, item.threshold.labels)
  }
  colnames(threshold.labels) <- sprintf("Item%d", 1:nrow(qmatrix))
  return(threshold.labels)
}

#' Obtain attribute profile probabilities from latent variable means
#'  
#' Gets all attribute profile probabilities given mu parameters (latent variable means)
#' @param mus a 1 by nclasses numeric vector of latent variable means
#' @return a numeric vector of all attribute profile probability values
GetClassProbsFromMus <- function(mus){
  nu <- sapply(mus, function(x) exp(x)/(sum(exp(mus))))
  names(nu) <- sprintf("nu_%d", 1:length(nu))
  return(nu)
}

#' Calculate value for pi.i.star.prime for one item
#' 
#' Given an NCRUM parameterization and values for pi.i.star, r.ia.star, and qmatrix, calculates the value for
#' for pi.i.star.prime for the given item
#' @param pi.i.star a numeric value for the probability that a person responds correctly given they have not mastered
#' any attributes before being penalized for not mastering an attribute.
#' @param r.ia.star a numeric 1 x nattributes vector of the penalty probabilities for this item of a correct response       
#' for each required attribute that has not been mastered.
#' @param q.ia a numeric 1 x nattributes vector of the qmatrix entries for the given item.
#' @return a numeric value representing the pi.i.star_prime for item i
CalculatePiIStarPrime <- function(pi.i.star, r.ia.star, q.ia){
  pi.i.star*prod(r.ia.star^(q.ia))
}

#' Calculate item probability of correct response for NCRUM parameterization
#'
#' Calculates the probability of a correct response for an item given an NCRUM parameterization
#' @param pi.i.star.prime a numeric value for the probability of a correct response for item i given 
#' that no attributes have been mastered.
#' @param r.ia.star a numeric 1 x nattributes vector of the penalty probabilities for this item of a correct response       
#' for each required attribute that has not been mastered.
#' @param q.ia a numeric 1 x nattributes vector of the qmatrix entries for the given item.
#' @param alphas a numeric 1 x nattributes vector with a 1 if a student has mastered an attribute and 0 otherwise.
#' @return the numeric value representing the probability of a correct response for the ith item.
GetProbCorrectNcrum <- function(pi.i.star.prime, r.ia.star, q.ia, alphas){
  pi.i.star.prime*prod(1/r.ia.star^(q.ia*alphas))
}

#' Calculate all item probabilities of correct response for NCRUM parameterization
#' 
#' Calculates the probabilities of a correct response for all items given an NCRUM parameterization
#' @param pi.i.star.prime a numeric 1 x nitems vector of the probabilities of a correct response for each item given 
#' that no attributes have been mastered.
#' @param r.ia.star.qmatrix a numeric nitems x nattributes vector of the penalty probabilities for each item of        
#' a correct response for each required attribute that has not been mastered.
#' @param qmatrix a numeric nitems by nattributes matrix that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param alphas a numeric 1 x n_attributes vector with a 1 if a student has mastered an attribute and 0 otherwise.
#' @return probabilities of a correct response for all items given a person's attribute profile.
GetAllProbsCorrectNcrum <- function(pi.i.star.prime, r.ia.star.qmatrix, qmatrix, alphas){
  sapply(1:nrow(qmatrix), function(x) GetProbCorrectNcrum(pi.i.star.prime[x], r.ia.star.qmatrix[x, ],
                                                          qmatrix[x, ], alphas))
}
#' Calculate item threshold values for NCRUM parameterization
#' 
#' Calculates all item threshold values for an NCRUM parameterization given parameter estimates,
#' qmatrix, attribute profiles, and threshold label matrix.
#' @param estimates0 a numeric vector of parameter estimates in the form of pis and rs
#' @param qmatrix a numeric nitems by nattributes matrix that specifies which items are required for mastery of each attribute (i.e., latent variable)
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles
#' @param threshold.labels an nclasses by nitems character matrix with appropriate item threshold labels
#' @return a numeric vector of item threshold values
GetThresholdValuesKernelPiR <- function(estimates0, qmatrix, pmatrix, threshold.labels){   
  pi.i.stars <- unlist(estimates0[grep('^pi', names(estimates0))])
  r.ia.stars <- unlist(estimates0[grep('^r', names(estimates0))])
  r.ia.star.qmatrix <- t(qmatrix)
  r.ia.star.qmatrix[r.ia.star.qmatrix==1] <- r.ia.stars
  r.ia.star.qmatrix <- t(r.ia.star.qmatrix)
  pi.i.star.primes <- sapply(1:nrow(qmatrix), function(x) CalculatePiIStarPrime(pi.i.stars[x], r.ia.star.qmatrix[x, ], qmatrix[x, ]))
  probs <- t(apply(pmatrix, 1, function(x) GetAllProbsCorrectNcrum(pi.i.star.primes, r.ia.star.qmatrix, qmatrix, x)))
  taus <- Logit(probs)
  taus <- unique(as.vector(t(taus)))
  names(taus) <- unique(as.vector(t(threshold.labels)))
  return(taus)
}

#' Calculate item threshold values for kernel parameterization
#'
#' Calculates the item threshold values for a kernel parameterization of a DCM model using lambdas
#' @param lambda.equations an nclasses by nitems character matrix with tau values
#' represented as a function of lambda values
#' @param estimates0 a numeric vector of parameter estimates in the form of lambdas and gammas
#' @param threshold.labels an nclasses by nitems character matrix with appropriate threshold labels.
#' @return a list of a numeric vector of item threshold values and an nclasses by nitems numeric matrix
#' of item threshold values
GetThresholdValuesKernel <- function(lambda.equations, estimates0, threshold.labels){
  nitems <- ncol(lambda.equations)
  nclasses <- nrow(lambda.equations)
  lambda.values <- unlist(estimates0[grep('^l', names(estimates0))])
  threshold.values <- matrix(0, nclasses, nitems)
  for (l in 1:length(lambda.values)){
    for (i in 1:nitems){
      if (any(grepl(names(lambda.values)[l], lambda.equations[,i]))){
        threshold.values[grep(names(lambda.values)[l], lambda.equations[,i]), i] <- threshold.values[grep(names(lambda.values)[l], lambda.equations[,i]), i] + lambda.values[l]
      }
    }
  }
  taus <- unique(as.vector(t(threshold.values)))
  names(taus) <- unique(as.vector(t(threshold.labels)))
  return(list(taus = taus, threshold.values = threshold.values))
}

#' Calcualte attribute means from gamma parameters
#' 
#' Give gamma parameters and all attribute profiles, calculates the attribute means (i.e., mus)
#' @param estimates0 a numeric vector of parameter estimates in the form of lambdas and gammas
#' @param pmatrix a numeric nclasses by nattributes matrix of all possible attribute profiles 
GetMusFromGammas <- function(estimates0, pmatrix){
  gamma.values <- unlist(estimates0[grep('^g', names(estimates0))])
  all.mus <- NULL
  all.mus <- c(all.mus, -1*(sum(gamma.values)))
  for (i in 2:(nrow(pmatrix) - 1)){
    cur.mu <- all.mus[1] + sum(gamma.values[which(pmatrix[i, ] == 1)])
    if (length(which(pmatrix[i, ] == 1)) > 1){
      tmp <- which(pmatrix[i, ] == 1)
      cur.mu <- cur.mu + gamma.values[grep(sprintf("g_%d_%s", length(tmp), paste(tmp, collapse="")), names(gamma.values))]
    }
    all.mus <- c(all.mus, cur.mu)
  }
  all.mus <- c(all.mus, 0)
  names(all.mus) <- sprintf("mu%d", 1:length(all.mus))
  return(all.mus)
}

#' GetThresholdValues
#' 
#' @param threshold.labels an nclasses by nitems character matrix with appropriate threshold labels.
#' @param taus numeric vector of taus correctly named 
#' threshold.labels.cls: String vector of names of thresholds for given class
GetThresholdValues <- function(threshold.labels, taus){
  t(sapply(1:nrow(threshold.labels), function(x) GetClassThresholds(taus, threshold.labels[x,])))
}

# Obtains correct threshold values ordered by item.
# inputs: taus: numeric vector of taus correctly named
#         threshold.labels.cls: String vector of names of thresholds for given class
# returns: A 1 by nitem vector of threshold values ordered by item.
GetClassThresholds <- function(taus, thresholds.labels.class){
  class.thresholds <- taus[names(taus) %in% thresholds.labels.class]
  class.thresholds[order(as.numeric(regmatches(names(class.thresholds)
                                               , regexpr('[0-9]+', names(class.thresholds)))))]
}
