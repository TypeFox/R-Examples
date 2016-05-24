################################################################################
# # MantelHaenszel
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Theoretical Mantel-Haenszel statistic under IRT assumptions
# # Description: Functions for calculating Theoretical Mantel-Haenszel statistic under IRT assumptions
# #
# # Inputs: NULL
# #
# # Outputs: functions
# #
# # File history:
# #   20120414: Creation
# #   20140424: Documentation adjusted to work with roxygen2. References
# #   added to documentation
################################################################################




################################################################################
# # Function CrossedProbabilities: Crossed products and Odds-ratio under
# # dichotomous IRT models
################################################################################

#' Calculates the crossed probabilities associated with the numerator and denominator of the odds-ratio under dichotomous IRT models
#'
#' @param thetaValue            A numeric value or array for the theta (ability) value(s) for which the odds will be calculated
#' @param itemParameters        A list containing the "focal" and "reference" item parameters. Item parameters  are assumed to be on the same scale.
#' @param irtModel              A string stating the irtModel used. May be one of "1pl", "2pl", or "3pl".
#' @param logistic              A logical indicating whether the logistic or the normal metric should be used.
#'
#' @return out A list containing the crossed products for the 'num' the numerator, 'den' the denominator for the odds-ratio, and 'or' the odds-ratio
#'
#' @references Roussos, L., Schnipke, D. & Pashley, P. (1999). A generalized formula for the Mantel-Haenszel Differential Item Functioning parameter. Journal of educational and behavioral statistics, 24(3), 293--322.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
CrossedProbabilities <- function(thetaValue, itemParameters, logistic, irtModel = "3pl") {

  if (!(irtModel %in% c("1pl", "2pl", "3pl"))) {
    stop("irtModel must be '1pl', '2pl' or '3pl'")
  }
  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  nItems <- nrow(itemParameters[["focal"]])

  if (irtModel == "1pl") {
    itemParameters[["focal"]] <- cbind(rep(1, nItems), itemParameters[["focal"]], rep(0, nItems))
    itemParameters[["reference"]] <- cbind(rep(1, nItems), itemParameters[["reference"]], rep(0, nItems))
  } else if (irtModel == "2pl") {
    itemParameters[["focal"]] <- cbind(itemParameters[["focal"]], rep(0, nItems))
    itemParameters[["reference"]] <- cbind(itemParameters[["reference"]], rep(0, nItems))
  }

  if (any(itemParameters[["focal"]][, 1] <= 0)) {
    stop("When discrimination parameters are included (2pl or 3pl), they must be in the first column of itemParameters and must be all positive")
  }

  if (any((itemParameters[["focal"]][, 3] < 0) || (itemParameters[["focal"]][, 3] > 1))) {
    stop("When guessing parameters are included (3pl), they must be in the third column of itemParameters and must be all in the interval [0, 1]")
  }
  if (any(itemParameters[["reference"]][, 1] <= 0)) {
    stop("When discrimination parameters are included (2pl or 3pl), they must be in the first column of itemParameters and must be all positive")
  }

  if (any((itemParameters[["reference"]][, 3] < 0) || (itemParameters[["reference"]][, 3] > 1))) {
    stop("When guessing parameters are included (3pl), they must be in the third column of itemParameters and must be all in the interval [0, 1]")
  }

  numProbabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters[["focal"]]))
  denProbabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters[["focal"]]))

  for (ii in seq(nrow(numProbabilities))) {
    numProbabilities[ii, ] <- (itemParameters[["reference"]][, 3] + ((1L - itemParameters[["reference"]][, 3]) *
                                                                     plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[["reference"]][, 2])),
                                                                            scale = 1L / (kD * itemParameters[["reference"]][, 1])))) *
(1L - (itemParameters[["focal"]][, 3] + ((1 - itemParameters[["focal"]][, 3]) *
                                         plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[["focal"]][, 2])),
                                                scale = 1L / (kD * itemParameters[["focal"]][, 1])))))

denProbabilities[ii, ] <-  (1L - (itemParameters[["reference"]][, 3] + ((1L - itemParameters[["reference"]][, 3]) *
                                                                        plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[["reference"]][, 2])),
                                                                               scale = 1L / (kD * itemParameters[["reference"]][, 1]))))) *
(itemParameters[["focal"]][, 3] + ((1 - itemParameters[["focal"]][, 3]) *
                                   plogis(q = thetaValue[ii], location = as.numeric(t(itemParameters[["focal"]][, 2])),
                                          scale = 1L / (kD * itemParameters[["focal"]][, 1]))))

  }

  out <- list(num = numProbabilities, den = denProbabilities, or = numProbabilities / denProbabilities)

  return(out)
}








################################################################################
# # Function IrtMh: Theoretical Mantel-Haenszel under dichotomous IRT
# # models
################################################################################

#' Calculates the Mantel-Haenszel theoretical parameter when a dichotomous IRT model holds
#'
#' @param itemParameters        A list containing the "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale.
#' @param irtModel              A string stating the irtModel used. May be one of "1pl", "2pl", or "3pl".
#' @param focalDistribution     A string stating the distribution assumed for the focal group.
#' @param focalDistrExtra       A list of extra parameters for the focal distribution function.
#' @param referenceDistribution A string stating the distribution assumed for the reference group.
#' @param referenceDistrExtra   A list of extra parameters for the reference distribution function.
#' @param groupRatio            A positive value indicating how many members of the reference group are expected for each member of the focal group.
#' @param logistic              A logical indicating whether the logistic or the normal metric should be used.
#' @param subdivisions          A numeric value stating the maximum number of subdivisions for adaptive quadrature.
#'
#' @return mh A list containing the asymptotic matrices for each item
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' threePlMh <- IrtMh(itemParameters = dichotomousItemParameters,  irtModel = "3pl",
#'                    focalDistribution = "norm", referenceDistribution = "norm",
#'                    focalDistrExtra = list(mean = 0, sd = 1),
#'                    referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1,
#'                    logistic = FALSE)
#'
#' @references Roussos, L., Schnipke, D. & Pashley, P. (1999). A generalized formula for the Mantel-Haenszel Differential Item Functioning parameter. Journal of educational and behavioral statistics, 24(3), 293--322.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
IrtMh <- function (itemParameters, irtModel = "2pl", focalDistribution = "norm", referenceDistribution = "norm",
                   focalDistrExtra = list(mean = 0, sd = 1), referenceDistrExtra = list(mean = 0, sd = 1),
                   groupRatio = 1, logistic = TRUE, subdivisions = 5000) {

  ################################################################################
  # # Auxiliary functions
  ################################################################################


  FocalDensity <- function (x, focalDistribution, focalDistrExtra) {
    pars <- focalDistrExtra
    pars$x <- x
    out  <- do.call(paste("d", focalDistribution, sep = ""), pars)
    return(out)
  }

  ReferenceDensity <- function (x, referenceDistribution, referenceDistrExtra) {
    pars <- referenceDistrExtra
    pars$x <- x
    out <- do.call(paste("d", referenceDistribution, sep = ""), pars)
    return(out)
  }

  DensityFactor <- function(x, focalDistribution, focalDistrExtra, referenceDistribution, referenceDistrExtra,
                            focalBaseRate, referenceBaseRate) {
    out <- (FocalDensity(x, focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra) *
            ReferenceDensity(x, referenceDistribution = referenceDistribution,
                             referenceDistrExtra = referenceDistrExtra))

    den <- ((focalBaseRate * FocalDensity(x, focalDistribution = focalDistribution,
                                          focalDistrExtra = focalDistrExtra)) +
            (referenceBaseRate * ReferenceDensity(x, referenceDistribution = referenceDistribution,
                                                  referenceDistrExtra = referenceDistrExtra)))

    out[out != 0] <- out[out != 0] / den[out != 0]

    return(out)
  }

  NumIntegrand <- function(x, itemParameters, focalBaseRate, referenceBaseRate, logistic, irtModel,
                           focalDistribution, focalDistrExtra,
                           referenceDistribution, referenceDistrExtra) {
    numCrossedProbabilities <- CrossedProbabilities(thetaValue = x, itemParameters = itemParameters,
                                                    logistic = logistic, irtModel = irtModel)$num

    densityFactor <- DensityFactor(x, focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra,
                                   referenceDistribution = referenceDistribution,
                                   referenceDistrExtra = referenceDistrExtra,
                                   focalBaseRate = focalBaseRate, referenceBaseRate = referenceBaseRate)

    out <- as.numeric(numCrossedProbabilities * densityFactor)

    return(out)
  }

  DenIntegrand <- function(x, itemParameters, focalBaseRate, referenceBaseRate, logistic, irtModel,
                           focalDistribution, focalDistrExtra,
                           referenceDistribution, referenceDistrExtra) {
    denCrossedProbabilities <- CrossedProbabilities(thetaValue = x, itemParameters = itemParameters,
                                                    logistic = logistic, irtModel = irtModel)$den

    densityFactor <- DensityFactor(x, focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra,
                                   referenceDistribution = referenceDistribution,
                                   referenceDistrExtra = referenceDistrExtra,
                                   focalBaseRate = focalBaseRate, referenceBaseRate = referenceBaseRate)

    out <- as.numeric(denCrossedProbabilities * densityFactor)

    return(out)
  }


  ################################################################################
  # # Base rate calculations
  ################################################################################

  focalBaseRate     <- 1 / (1 + groupRatio)
  referenceBaseRate <- groupRatio / (1 + groupRatio)
  nItems            <- nrow(itemParameters[["focal"]])

  num <- numeric(nItems)
  den <- numeric(nItems)

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }
  if (irtModel == "1pl") {
    mh <- as.numeric(exp(-kD * (itemParameters[["reference"]] - itemParameters[["focal"]])))
  } else {
    for (ii in seq(nItems)) {
      iiItemParameters <- list(focal = matrix(itemParameters[["focal"]][ii, ], nrow = 1),
                               reference = matrix(itemParameters[["reference"]][ii, ], nrow = 1))
      if (irtModel == "2pl" & (iiItemParameters[["focal"]][, 1] == iiItemParameters[["reference"]][, 1])) {
        num[ii] <- exp(-kD * iiItemParameters[["focal"]][, 1] *
                       (iiItemParameters[["reference"]][, 2] - iiItemParameters[["focal"]][, 2]))
        den[ii] <- 1
      } else {
        num[ii] <- integrate(f = NumIntegrand, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                             focalBaseRate, referenceBaseRate, irtModel = irtModel,
                             itemParameters = iiItemParameters, logistic = logistic,
                             focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra,
                             referenceDistribution = referenceDistribution, referenceDistrExtra = referenceDistrExtra)$value
        den[ii] <- integrate(f = DenIntegrand, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                             focalBaseRate, referenceBaseRate, irtModel = irtModel,
                             itemParameters = iiItemParameters, logistic = logistic,
                             focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra,
                             referenceDistribution = referenceDistribution,
                             referenceDistrExtra = referenceDistrExtra)$value
      }
    }
    mh <- num / den
  }

  return(mh)
}








################################################################################
# # Function DeltaMhIrt: Tranform the MH statistic to the ETS Delta
# # metric
################################################################################

#' Obtains the ETS Delta measure for Mantel-Haneszel DIF statistic effect size.
#'
#' @param mh A numeric vector containing the MH statistic values
#' @param logistic              A logical indicating whether the logistic or the normal metric should be used.
#'
#' @return delta A numeric vector containing the delta values
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' threePlMh <- IrtMh(itemParameters = dichotomousItemParameters,  irtModel = "3pl",
#'                    focalDistribution = "norm",
#'                    referenceDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
#'                    referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1,
#'                    logistic = FALSE)
#' delta3pl <- DeltaMhIrt(threePlMh)
#'
#' @references Holland, P.W., and Thayer, D.T. (1988). Differential Item Performance and the Mantel-Haenszel Procedure. In H. Wainer and H.I. Braun (Eds.), Test Validity. Hillsdale, NJ: Erlbaum.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
DeltaMhIrt <- function (mh, logistic = FALSE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  } 

  delta <- (-4 / kD) * log(mh)

  return(delta)
}
