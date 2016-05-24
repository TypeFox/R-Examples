################################################################################
# # IRTSE.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Analytical Standard Errors for dichotomous and polytomous IRT Models
# # Description: Functions for calculating analytical standard errors for dichotomous and polytomous IRT models
# #              including item parameters' estimation covariances
# #
# # Inputs: None
# #
# # Outputs: functions
# #
# # File history:
# #   20120425: Creation
# #   20140421: Documentation adjusted to work with roxygen2. References
# #   added to documentation
################################################################################




################################################################################
# # Function ProductProbabilities: Product of item response
# # probabilities for dichotomous IRT models
################################################################################

#' Calculates the product of item response probabilities for dichotomous IRT models
#'
#' @param thetaValue            A numeric value or array for the theta (ability) value(s) for which the product will be calculated
#' @param itemParameters        A matrix containing the item parameters.
#' @param irtModel              A string stating the irtModel used. May be one of "1pl", "2pl", or "3pl".
#' @param logistic              A logical indicating whether the logistic or the normal metric should be used.
#'
#' @return pq A numeric matrix containing the crossed product on each thetaValue for each item.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
ProductProbabilities <- function(thetaValue, itemParameters, logistic, irtModel = "3pl") {

  if (!(irtModel %in% c("1pl", "2pl", "3pl"))) {
    stop("irtModel must be '1pl', '2pl' or '3pl'")
  }

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  nItems <- nrow(itemParameters)

  if (irtModel == "1pl") {
    itemParameters <- cbind(rep(1, nItems), itemParameters, rep(0, nItems))
  } else if (irtModel == "2pl") {
    itemParameters <- cbind(itemParameters, rep(0, nItems))
  }

  if (any(itemParameters[, 1] <= 0)) {
    stop("When discrimination parameters are included (2pl or 3pl), they must be in the first column of itemParameters and must be all positive")
  }

  if (any((itemParameters[, 3] < 0) || (itemParameters[, 3] > 1))) {
    stop("When guessing parameters are included (3pl), they must be in the third column of itemParameters and must be all in the interval [0, 1]")
  }

  pq <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters))

  for (ii in seq(nrow(pq))) {
    pq[ii, ] <- ((itemParameters[, 3] + ((1L - itemParameters[, 3]) *
                                         plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[, 2])),
                                                scale = 1L / (kD * itemParameters[, 1])))) *
                 (1L - (itemParameters[, 3] + ((1 - itemParameters[, 3]) *
                                               plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[, 2])),
                                                      scale = 1L / (kD * itemParameters[, 1]))))))
  }

  return(pq)
}








################################################################################
# # Function AseIrt: Asymptotic covariance matrices of item parameter
# # estimates
################################################################################

#' Calculates the asymptotic covariance matrices for item parameters according with the IRT model.
#'
#' @param itemParameters         A matrix or vector containing the item difficulties.
#' @param distribution           A string character indicating the generic name for the assumed distribution. Defaults to 'norm' for normal distribution.
#' @param distributionParameters A list of extra parameters for the distribution function.
#' @param logistic               A logical indicating whether the logistic or the normal metric should be used.
#' @param sampleSize             A value indicating the sample size.
#' @param irtModel               A string stating the IRT model for all items.
#' @param subdivisions           A numeric value stating the maximum number of subdivisions for adaptive quadrature.
#'
#' @return ase A list containing the asymptotic matrices for each item
#'
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically derived standard errors of Item Response Theory item parameter estimates. Journal of educational measurement, 41(2), 85--117.
#'
#' @export
#'
#' @examples
#' # # Not run
#' # #
#' # # data(dichotomousItemParameters)
#' # # threePlAse <- list()
#' # # threePlAse[['focal']] <- AseIrt(itemParameters = dichotomousItemParameters[['focal']],
#' # #                                 logistic = TRUE, sampleSize = 500, irtModel = '3pl')
#' # # threePlAse[['reference']] <- AseIrt(itemParameters = dichotomousItemParameters[['reference']],
#' # #                                     logistic = TRUE, sampleSize = 500, irtModel = '3pl')
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
AseIrt <- function (itemParameters, distribution = "norm", distributionParameters = list(mean = 0, sd = 1),
                    logistic = TRUE, sampleSize = 1, irtModel = "3pl", subdivisions = 5000) {
  # # TODO: Change so that each item may be modelled by a different IRT
  # # model

  CalculateAse <- switch(irtModel,
                         "1pl" = Ase1pl,
                         "2pl" = Ase2pl,
                         "3pl" = Ase3pl,
                         "grm" = stop("Not yet implemented"),
                         "pcm" = stop("Not yet implemented"),
                         stop("irtModel not known or not implemented"))

  ase <- CalculateAse(itemParameters = itemParameters, distribution = distribution,
                      distributionParameters = distributionParameters, logistic = logistic,
                      sampleSize = sampleSize, subdivisions = subdivisions)

  return(ase)
}








################################################################################
# # Function Ase1pl:  Asymptotic variances of item parameter for the 1PL
# # IRT model
################################################################################

#' Calculates the asymptotic variance for difficulty parameter estimates under the 1pl model
#'
#' @param itemParameters         A matrix or vector containing the item difficulties.
#' @param distribution           A string character indicating the generic name for the assumed distribution.
#' @param distributionParameters A list of extra parameters for the distribution function.
#' @param logistic               A logical indicating whether the logistic or the normal metric should be used.
#' @param sampleSize             A value indicating the sample size.
#' @param subdivisions           A numeric value stating the maximum number of subdivisions for adaptive quadrature.
#'
#' @return ase A list containing the asymptotic variances for each item
#'
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically derived standard errors of Item Response Theory item parameter estimates. Journal of educational measurement, 41(2), 85--117.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Ase1pl <- function (itemParameters, distribution = "norm", distributionParameters = list(mean = 0, sd = 1),
                    logistic = TRUE, sampleSize = 1, subdivisions = 5000) {

  Integrand <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars         <- distributionParameters
    pars$x       <- x
    densityValue <- do.call(paste("d", distribution, sep = ""), pars)
    pq           <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "1pl")

    out <- pq * densityValue
    return(out)
  }

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  nItems <- length(itemParameters)
  ase    <- list()

  for (ii in 1:nItems) {
    iiItemParameters <- matrix(itemParameters[ii], nrow = 1)
    ase[[ii]] <- 1 / (sampleSize * kD * kD * integrate(f = Integrand, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                                       distribution = distribution, distributionParameters = distributionParameters,
                                                       itemParameters = iiItemParameters, logistic = logistic)$value)
  }

  return(ase)
}








################################################################################
# # Function Ase2pl:  Asymptotic covariances of item parameter for the 2PL
# # IRT model
################################################################################

#' Calculates the asymptotic covariance matrix of item parameter estimates under the 2pl model
#'
#' @param itemParameters         A matrix or vector containing the item difficulties.
#' @param distribution           A string character indicating the generic name for the assumed distribution.
#' @param distributionParameters A list of extra parameters for the distribution function.
#' @param logistic               A logical indicating whether the logistic or the normal metric should be used.
#' @param sampleSize             A value indicating the sample size.
#' @param subdivisions           A numeric value stating the maximum number of subdivisions for adaptive quadrature.
#'
#' @return ase A list containing the asymptotic matrices for each item
#'
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically derived standard errors of Item Response Theory item parameter estimates. Journal of educational measurement, 41(2), 85--117.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Ase2pl <- function (itemParameters, distribution = "norm", distributionParameters = list(mean = 0, sd = 1),
                    logistic = TRUE, sampleSize = 1, subdivisions = 5000) {

  IntegrandDiscrimination <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "2pl")
    thetaDifficuty <- (x - itemParameters[, 2])

    out <- thetaDifficuty * thetaDifficuty * pq * densityValue
    return(out)
  }

  IntegrandDiscriminationDifficulty <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "2pl")
    thetaDifficuty <- (x - itemParameters[, 2])

    out <- thetaDifficuty * pq * densityValue
    return(out)
  }

  IntegrandDifficulty <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars         <- distributionParameters
    pars$x       <- x
    densityValue <- do.call(paste("d", distribution, sep = ""), pars)
    pq           <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "2pl")

    out <- pq * densityValue
    return(out)
  }

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  nItems <- nrow(itemParameters)
  ase    <- list()

  for (ii in 1:nItems) {
    iiItemParameters <- matrix(itemParameters[ii, ], nrow = 1)
    ase[[ii]] <- matrix(nrow = 2, ncol = 2)
    ase[[ii]][1, 1] <- (sampleSize * kD * kD *
                        integrate(f = IntegrandDiscrimination, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                  distribution = distribution, distributionParameters = distributionParameters,
                                  itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][1, 2] <- ase[[ii]][2, 1]  <- (-sampleSize * kD * kD * iiItemParameters[, 1] *
                                            integrate(f = IntegrandDiscriminationDifficulty, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                                      distribution = distribution, distributionParameters = distributionParameters,
                                                      itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][2, 2] <- (sampleSize * kD * kD * iiItemParameters[, 1] * iiItemParameters[, 1] *
                        integrate(f = IntegrandDifficulty, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                  distribution = distribution, distributionParameters = distributionParameters,
                                  itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]] <- solve(ase[[ii]])
  }

  return(ase)
}








################################################################################
# # Function Ase3pl:  Asymptotic covariances of item parameter for the 3PL
# # IRT model
################################################################################

#' Calculates the asymptotic covariance matrix of item parameter estimates under the 3pl model
#'
#' @param itemParameters         A matrix or vector containing the item difficulties.
#' @param distribution           A string character indicating the generic name for the assumed distribution.
#' @param distributionParameters A list of extra parameters for the distribution function.
#' @param logistic               A logical indicating whether the logistic or the normal metric should be used.
#' @param sampleSize             A value indicating the sample size.
#' @param subdivisions           A numeric value stating the maximum number of subdivisions for adaptive quadrature.
#'
#' @return ase A list containing the asymptotic matrices for each item
#'
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically derived standard errors of Item Response Theory item parameter estimates. Journal of educational measurement, 41(2), 85--117.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Ase3pl <- function (itemParameters, distribution = "norm", distributionParameters = list(mean = 0, sd = 1),
                    logistic = TRUE, sampleSize = 1, subdivisions = 5000) {

  kUsualGuessing <- 0.2 # #

  IntegrandDiscrimination <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pqAst          <- ProductProbabilities(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic, irtModel = "2pl")
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")
    thetaDifficuty <- (x - itemParameters[, 2])

    out <- thetaDifficuty * thetaDifficuty * pqAst * pqAst * densityValue
    out[pq != 0] <- out[pq != 0] / pq[pq != 0]
    out[!is.finite(out)] <- 0
    return(out)
  }

  IntegrandDiscriminationDifficulty <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pqAst          <- ProductProbabilities(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic, irtModel = "2pl")
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")
    thetaDifficuty <- (x - itemParameters[, 2])

    out <- thetaDifficuty * pqAst * pqAst * densityValue
    out[out != 0] <- out[out != 0] / pq[out != 0]
    out[!is.finite(out)] <- 0
    return(out)
  }

  IntegrandDiscriminationGuessing <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pqAst          <- ProductProbabilities(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic, irtModel = "2pl")
    qAst           <- 1 - Calculate2plProb(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic)
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")
    thetaDifficuty <- (x - itemParameters[, 2])

    out <- thetaDifficuty * pqAst * qAst * densityValue
    out[out != 0] <- out[out != 0] / pq[out != 0]
    out[!is.finite(out)] <- 0
    return(out)
  }

  IntegrandDifficulty <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars         <- distributionParameters
    pars$x       <- x
    densityValue <- do.call(paste("d", distribution, sep = ""), pars)
    pqAst          <- ProductProbabilities(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic, irtModel = "2pl")
    pq           <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")

    out <- pqAst * pqAst * densityValue
    out[out != 0] <- out[out != 0] / pq[out != 0]
    out[!is.finite(out)] <- 0
    return(out)
  }

  IntegrandDifficultyGuessing <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    pqAst          <- ProductProbabilities(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic, irtModel = "2pl")
    qAst           <- 1 - Calculate2plProb(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic)
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")

    out <- pqAst * qAst * densityValue
    out[out != 0] <- out[out != 0] / pq[out != 0]
    out[!is.finite(out)] <- 0
    return(out)
  }

  IntegrandGuessing <- function (x, itemParameters, distribution, distributionParameters, logistic) {
    pars           <- distributionParameters
    pars$x         <- x
    densityValue   <- do.call(paste("d", distribution, sep = ""), pars)
    qAst           <- 1 - Calculate2plProb(thetaValue = x, itemParameters = matrix(itemParameters[, -3], ncol = 2), logistic = logistic)
    pq             <- ProductProbabilities(thetaValue = x, itemParameters = itemParameters, logistic = logistic, irtModel = "3pl")

    out <- qAst * qAst * densityValue
    out[out != 0] <- out[out != 0] / pq[out != 0]
    out[!is.finite(out)] <- 1
    return(out)
  }

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  nItems <- nrow(itemParameters)
  ase    <- list()

  for (ii in 1:nItems) {
    iiItemParameters <- matrix(itemParameters[ii, ], nrow = 1)
    ase[[ii]] <- matrix(nrow = 3, ncol = 3)
    ase[[ii]][1, 1] <- (sampleSize * kD * kD * (1 - iiItemParameters[, 3]) * (1 - iiItemParameters[, 3]) *
                        integrate(f = IntegrandDiscrimination, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                  distribution = distribution, distributionParameters = distributionParameters,
                                  itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][1, 2] <- ase[[ii]][2, 1]  <- (-sampleSize * kD * kD * iiItemParameters[, 1] * (1 - iiItemParameters[, 3]) * (1 - iiItemParameters[, 3]) *
                                            integrate(f = IntegrandDiscriminationDifficulty, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                                      distribution = distribution, distributionParameters = distributionParameters,
                                                      itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][1, 3] <- ase[[ii]][3, 1]  <- (sampleSize * kD * (1 - iiItemParameters[, 3]) *
                                            integrate(f = IntegrandDiscriminationGuessing, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                                      distribution = distribution, distributionParameters = distributionParameters,
                                                      itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][2, 2] <- (sampleSize * kD * kD * iiItemParameters[, 1] * iiItemParameters[, 1] * (1 - iiItemParameters[, 3]) * (1 - iiItemParameters[, 3]) *
                        integrate(f = IntegrandDifficulty, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                  distribution = distribution, distributionParameters = distributionParameters,
                                  itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][2, 3] <- ase[[ii]][3, 2]  <- (-sampleSize * kD * iiItemParameters[, 1] * (1 - iiItemParameters[, 3]) *
                                            integrate(f = IntegrandDifficultyGuessing, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                                      distribution = distribution, distributionParameters = distributionParameters,
                                                      itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]][3, 3] <- (sampleSize *
                        integrate(f = IntegrandGuessing, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                  distribution = distribution, distributionParameters = distributionParameters,
                                  itemParameters = iiItemParameters, logistic = logistic)$value)

    ase[[ii]] <- solve(ase[[ii]])

    if (qnorm(0.975) * sqrt(ase[[ii]][3, 3]) >= kUsualGuessing) {
      warning("Sample size may be to small for stable guessing estimates for item ", ii)
    }
  }

  return(ase)
}
