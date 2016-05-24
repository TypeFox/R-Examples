################################################################################
# # NCDIF.R
# # R Versions: 2.14.1
# #
# # Author(s): Victor H. Cervantes
# #
# # General NCDIF function
# # Description: Functions for calculating NCDIF
# #
# # Inputs: NULL
# #
# # Outputs: functions
# #
# # File history:
# #   20120304: Creation
# #   20140421: Documentation adjusted to work with roxygen2. References
# #   added to documentation
# #   20140518: Updated PlotNcdif for latest ggplot2 compability
# #   20140518: Changed ... option for indices functions
################################################################################




################################################################################
# # Function Ncdif: NonCompensatory Differential Item Functioning index
################################################################################

#' Calculates NCDIF index for an item with given item parameters of focal and reference groups.
#'
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalAbilities    If NULL, NCDIF is calculated by numerical integration of focal distribution. If not NULL, it must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution A string stating the distribution name to be used for integrating. Only used if focalAbilities is NULL.
#' @param focalDistrExtra   Extra parameters for the focal group distribution function if needed.
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities is NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return ncdif Numeric vector with the NCDIF index value for each item.
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' threePlNcdif <- Ncdif(itemParameters = dichotomousItemParameters, irtModel = '3pl',
#'                       focalAbilities = NULL, focalDistribution = "norm",
#'                       subdivisions = 5000, logistic = TRUE)
#'
#' @references Raju, N. S., van der Linden, W. J., & Fleer, P. F. (1995). An IRT-based internal measure of test bias with applications for differential item functioning. Applied Psychological Measurement, 19, 353--368.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Ncdif <- function (itemParameters, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm",
                   subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(mean = 0, sd = 1)) {

  NcdifPoint <- function (x, itemParameters, focalDistribution, irtModel, logistic, focalDistrExtra) {
    itemDifference <- CalculateItemDifferences(thetaValue = x, itemParameters = itemParameters, irtModel = irtModel,
                                               logistic = logistic)
    pars <- focalDistrExtra
    pars$x <- x
    focalGroupDensity <- do.call(paste("d", focalDistribution, sep = ""), pars)

    out <- focalGroupDensity * ((itemDifference) ^ 2)
    return(out)
  }

  if (is.null(focalAbilities)) {
    nItems <- nrow(itemParameters[["focal"]])
    ncdif  <- numeric(nItems)

    for (ii in 1:nItems) {
      ncdif[ii] <- integrate(f = NcdifPoint, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                             itemParameters = list(focal     = matrix(itemParameters[["focal"]][ii, ], nrow = 1),
                                                   reference = matrix(itemParameters[["reference"]][ii, ], nrow = 1)),
                             focalDistribution = focalDistribution,
                             focalDistrExtra = focalDistrExtra, irtModel = irtModel,
                             logistic = logistic)$value
    }
  } else {
    ncdif <- colMeans((CalculateItemDifferences(thetaValue = focalAbilities, itemParameters = itemParameters,
                                                irtModel = irtModel, logistic = logistic) ^ 2))
  }

  ncdif <- as.numeric(ncdif)

  return(ncdif)
}








################################################################################
# #  Function Cdif: Compensatory Differential Item Functioning index
################################################################################

#' Calculates CDIF index for an item with given item parameters of focal and reference groups.
#'
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalAbilities    If NULL, NCDIF is calculated by numerical integration of focal distribution. If not NULL, it must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution A string stating the distribution name to be used for integrating. Only used if focalAbilities is NULL.
#' @param focalDistrExtra   Extra parameters for the focal group distribution function if needed.
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities is NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return cdif Numeric vector with the CDIF index value for each item.
#'
#' @export
#'
#' @examples
#'
#' # # Not run
#'
#' # # data(dichotomousItemParameters)
#' # # threePlCdif <- Cdif(itemParameters = dichotomousItemParameters, irtModel = '3pl',
#' # #                     focalAbilities = NULL, focalDistribution = "norm",
#' # #                     subdivisions = 5000, logistic = TRUE)
#'
#' @references Raju, N. S., van der Linden, W. J., & Fleer, P. F. (1995). An IRT-based internal measure of test bias with applications for differential item functioning. Applied Psychological Measurement, 19, 353--368.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Cdif <- function (itemParameters, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm",
                  subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(mean = 0, sd = 1)) {

  CdifPoint <- function (x, itemParameters, iiItem, jjItem, focalDistribution, irtModel, logistic,
                         focalDistrExtra) {
    iiItemDifference <- CalculateItemDifferences(thetaValue = x,
                                                 itemParameters = list(focal     = matrix(itemParameters[["focal"]][iiItem, ],     nrow = 1),
                                                                       reference = matrix(itemParameters[["reference"]][iiItem, ], nrow = 1)),
                                                 irtModel = irtModel,
                                                 logistic = logistic)
    jjItemDifference <- CalculateItemDifferences(thetaValue = x,
                                                 itemParameters = list(focal     = matrix(itemParameters[["focal"]][jjItem, ],     nrow = 1),
                                                                       reference = matrix(itemParameters[["reference"]][jjItem, ], nrow = 1)),
                                                 irtModel = irtModel,
                                                 logistic = logistic)

    pars <- focalDistrExtra
    pars$x <- x
    focalGroupDensity <- do.call(paste("d", focalDistribution, sep = ""), pars)

    out <- focalGroupDensity * iiItemDifference * jjItemDifference

    return(out)
  }

  if (is.null(focalAbilities)) {
    nItems <- nrow(itemParameters[["focal"]])
    cdif   <- matrix(nrow = nItems, ncol = nItems)

    for (ii in 1:nItems) {
      for (jj in 1:nItems) {
        if (ii <= jj) {
          cdif[ii, jj] <- integrate(f = CdifPoint, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                                    itemParameters = itemParameters,
                                    iiItem = ii, jjItem = jj,
                                    focalDistribution = focalDistribution,
                                    focalDistrExtra = focalDistrExtra, irtModel = irtModel,
                                    logistic = logistic)$value
        }
      }
    }

    cdif[lower.tri(cdif)] <- t(cdif)[lower.tri(cdif)]
    cdif <- rowSums(cdif)
  } else {
    nIndividuals <- length(focalAbilities)
    cdif <- CalculateItemDifferences(thetaValue = focalAbilities, itemParameters = itemParameters,
                                     irtModel = irtModel, logistic = logistic)
    cdif <- rowSums((t(cdif) %*% cdif) / nIndividuals)
  }

  cdif <- as.numeric(cdif)

  return(cdif)
}








################################################################################
# #  Function Dtf: Differential Test Functioning index
################################################################################

#' Calculates DTF index for a set of items with given item parameters of focal and reference groups.
#'
#' @param cdif              A numeric vector of CDIF values for the test items. If NULL it is calculated using itemParameters and the other arguments.
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Only used if cdif is NULL. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm". Only used if cdif is NULL.
#' @param focalAbilities    If NULL, CDIF is calculated by numerical integration of focal distribution. If not NULL, it must be a numerical vector containing the abilities for the individuals in the focal group. Only used if cdif is NULL.
#' @param focalDistribution A string stating the distribution name to be used for integrating. Only used if focalAbilities and cdif are NULL.
#' @param focalDistrExtra   Extra parameters for the focal group distribution function if needed.
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities and cdif are NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used. Only used if cdif is NULL.
#'
#' @return dtf Numeric vector with the CDIF index value for each item.
#'
#' @export
#'
#' @examples
#'
#' # # Not run
#'
#' # # data(dichotomousItemParameters)
#' # # threePlCdif <- Cdif(itemParameters = dichotomousItemParameters, irtModel = '3pl',
#' # #                     focalAbilities = NULL, focalDistribution = "norm",
#' # #                     subdivisions = 5000, logistic = TRUE)
#' # # threePlDtf  <- Dtf(cdif = threePlCdif)
#'
#' @references Raju, N. S., van der Linden, W. J., & Fleer, P. F. (1995). An IRT-based internal measure of test bias with applications for differential item functioning. Applied Psychological Measurement, 19, 353--368.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Dtf <- function (cdif = NULL, itemParameters = NULL, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm",
                 subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(mean = 0, sd = 1)) {

  if (is.null(cdif) & is.null(itemParameters)) {
    stop("Either the CDIF values for each item or the item parameters for calculating them must be supplied")
  }

  if (is.null(cdif)) {
    cdif <- Cdif(itemParameters = itemParameters, irtModel = irtModel, focalAbilities = focalAbilities,
                 focalDistribution = focalDistribution, subdivisions = subdivisions, logistic = logistic,
                 focalDistrExtra)
  }

  dtf <- sum(cdif)

  return(dtf)
}








################################################################################
# #  Function CalculateItemDifferences: Calculate differences between
# #  two item option characteristic curves
################################################################################

#' Calculates the differences between two item option characteristic curves for all options (minus one).
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) for which the difference will be calculated
#' @param itemParameters A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel       A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return difference A numeric matrix with the differences on probabilities or on expected score for each item between focal and reference groups.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
CalculateItemDifferences <- function (thetaValue, itemParameters, irtModel = "2pl", logistic = TRUE) {

  # # Todo: Change to allow each item to be modeled by a different IRT
  # # model

  if (!is.numeric(thetaValue)) {
    stop("thetaValue must be a numeric value or array")
  }

  CalculateProb <- switch(irtModel,
                          "1pl" = Calculate1plProb,
                          "2pl" = Calculate2plProb,
                          "3pl" = Calculate3plProb,
                          "grm" = CalculateGrmExp,
                          "pcm" = CalculatePcmExp,
                          stop("irtModel not known or not implemented"))

  focalProb     <- CalculateProb(thetaValue = thetaValue,
                                 itemParameters = itemParameters[["focal"]],
                                 logistic = logistic)

  referenceProb <- CalculateProb(thetaValue = thetaValue,
                                 itemParameters = itemParameters[["reference"]],
                                 logistic = logistic)


  difference <- focalProb - referenceProb

  return(difference)
}








################################################################################
# #  Function Calculate1plProb: Calculate the item success probability under the 1PL model.
################################################################################

#' Calculates the item success probability under the 1PL model.
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) where the difference will be calculated
#' @param itemParameters A vector or column matrix containing the numeric values of item difficulties
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return probabilities A numeric matrix with the probabilities on each thetaValue for each item.
#'
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Calculate1plProb <- function (thetaValue, itemParameters, logistic = TRUE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  if (length(dim(itemParameters)) > 0) {
    itemParameters <- as.numeric(itemParameters)
  }

  probabilities <- matrix(nrow = length(thetaValue), ncol = length(itemParameters))

  for (ii in seq(nrow(probabilities))) {
    probabilities[ii, ] <- plogis(q = thetaValue[ii],
                                  location = itemParameters,
                                  scale = 1 / kD)
  }

  return(probabilities)
}








################################################################################
# #  Function Calculate2plProb: Calculate the item success probability under the 2PL model.
################################################################################

#' Calculates the item success probability under the 2PL model.
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) where the difference will be calculated
#' @param itemParameters A matrix containing the numeric values of item discriminations on the first column and item difficulties on the second
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return probabilities A numeric matrix with the probabilities on each thetaValue for each item.
#'
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Calculate2plProb <- function (thetaValue, itemParameters, logistic = TRUE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  if (any(itemParameters[, 1] <= 0)) {
    stop("Discrimination parameters must be in the first column of itemParameters and must be all positive")
  }

  probabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters))

  for (ii in seq(nrow(probabilities))) {
    probabilities[ii, ] <- plogis(q = thetaValue[ii],
                                  location = as.numeric(t(itemParameters[, 2])),
                                  scale = 1 / (kD * itemParameters[, 1]))
  }

  return(probabilities)
}








################################################################################
# #  Function Calculate3plProb: Calculate the item success probability under the 3PL model.
################################################################################

#' Calculates the item success probability under the 3PL model.
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) where the difference will be calculated
#' @param itemParameters A matrix containing the numeric values of item discriminations on the first column, item difficulties on the second and item guessing parameters on the third
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return probabilities A numeric matrix with the probabilities on each thetaValue for each item.
#'
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Calculate3plProb <- function (thetaValue, itemParameters, logistic = TRUE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  if (any(itemParameters[, 1] <= 0)) {
    stop("Discrimination parameters must be in the first column of itemParameters and must be all positive")
  }

  if (any((itemParameters[, 3] < 0) || (itemParameters[, 3] > 1))) {
    stop("Guessing parameters must be in the third column of itemParameters and must be all in the interval [0, 1]")
  }

  if (is.data.frame(itemParameters)) {
    itemParameters <- as.matrix(itemParameters)
  }

  probabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters))

  for (ii in seq(nrow(probabilities))) {
    probabilities[ii, ] <-  itemParameters[, 3] + ((1 - itemParameters[, 3]) * plogis(q = thetaValue[ii],
                                                                                      location = as.numeric(t(itemParameters[, 2])),
                                                                                      scale = 1 / (kD * itemParameters[, 1])))
  }

  return(probabilities)
}








################################################################################
# #  Function CalculateGrmExp: Calculate the expected item score under the GRM model.
################################################################################

#' Calculates the expected item score under the GRM model.
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) where the difference will be calculated
#' @param itemParameters A matrix containing the numeric values of item discriminations on the first column and category thresholds on the rest columns where the (column position - 1) indicates the category score or weight.
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return expectedScore A numeric matrix with the expected score on each thetaValue for each item.
#'
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#' @references Oshima, T. & Morris, S. (2008). Raju's Differential Functioning of Items and Tests (DFIT). Educational Measurement: Issues and Practice, 27(3), 43--50.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
CalculateGrmExp <- function (thetaValue, itemParameters, logistic = TRUE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  if (any(itemParameters[, 1] <= 0)) {
    stop("Discrimination parameters must be in the first column of itemParameters and must be all positive")
  }

  if (is.data.frame(itemParameters)) {
    itemParameters <- as.matrix(itemParameters)
  }

  probabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters) * (ncol(itemParameters) - 1))
  itemIndices   <- matrix(rep(1:nrow(itemParameters), each = length(thetaValue) * (ncol(itemParameters) - 1)),
                          nrow = length(thetaValue), ncol = nrow(itemParameters) * (ncol(itemParameters) - 1))

  for (ii in seq(nrow(probabilities))) {
    probabilities[ii, ] <- plogis(q = thetaValue[ii],
                                  location = as.numeric(t(itemParameters[, -1])),
                                  scale = 1 / (kD * rep(itemParameters[, 1], each = ncol(itemParameters) - 1)))
  }

  expectedScore <- tapply(probabilities, itemIndices,
                          function (y)
                            apply(matrix(y, nrow(probabilities)), 1,
                                  function (x)
                                    sum((length(x) + 1):1 * diff(c(0, x[length(x):1], 1)), na.rm = TRUE)))

  expectedScore <- as.matrix(as.data.frame.list(expectedScore))

  return(expectedScore)
}








################################################################################
# #  Function CalculatePcmExp: Calculate the expected item score under the (G)PCM model.
################################################################################

#' Calculates the expected item score under the (G)PCM model.
#'
#' @param thetaValue     A numeric value or array for the theta (ability) value(s) where the difference will be calculated
#' @param itemParameters A matrix containing the numeric values of item discriminations on the first column and category thresholds on the rest columns where the (column position - 1) indicates the category score or weight.
#' @param logistic       A logical value stating if the IRT model will use the logistic or the normal metric.
#'
#' @return expectedScore A numeric matrix with the expected score on each thetaValue for each item.
#'
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#' @references Oshima, T. & Morris, S. (2008). Raju's Differential Functioning of Items and Tests (DFIT). Educational Measurement: Issues and Practice, 27(3), 43--50.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
CalculatePcmExp <- function (thetaValue, itemParameters, logistic = TRUE) {

  if (logistic) {
    kD <- 1
  } else {
    kD <- 1.702
  }

  if (any(itemParameters[, 1] <= 0)) {
    stop("Discrimination parameters must be in the first column of itemParameters and must be all positive")
  }

  if (is.data.frame(itemParameters)) {
    itemParameters <- as.matrix(itemParameters)
  }

  probabilities <- matrix(nrow = length(thetaValue), ncol = nrow(itemParameters) * (ncol(itemParameters) - 1))
  itemIndices   <- matrix(rep(1:nrow(itemParameters), each = length(thetaValue) * (ncol(itemParameters) - 1)),
                          nrow = length(thetaValue), ncol = nrow(itemParameters) * (ncol(itemParameters) - 1))

  for (ii in seq(nrow(probabilities))) {
    probabilities[ii, ] <- kD * rep(itemParameters[, 1], each = ncol(itemParameters) - 1) *
                           (thetaValue[ii] - as.numeric(t(itemParameters[, -1])))
  }

  expectedScore <- tapply(probabilities, itemIndices,
                          function (y)
                            apply(matrix(y, nrow(probabilities)), 1,
                                  function (x)
                                    sum(1:(length(x) + 1) * exp(c(0, cumsum(x))) / sum(exp(c(0, cumsum(x)))))
                                  ))

  expectedScore <- as.matrix(as.data.frame.list(expectedScore))

  isNanExp <- is.nan(expectedScore)
  isInfExp <- expectedScore == Inf
  expectedScore[isNanExp | isInfExp] <- apply(itemParameters, 1, function (x) sum(!is.na(x)))

  return(expectedScore)
}








################################################################################
# #  Function PlotNcdif: Plot the item characteristic (expected score)
# #  curve for focal and reference groups for the iiItem along with a
# #  representation of the focal group density.
################################################################################
if (getRversion() >= "2.15.1") utils::globalVariables(c("icc", "expected"))

#' Plot the item characteristic (expected score) curve for focal and reference groups for the iiItem along with a
#'   representation of the focal group density.
#'
#' @param iiItem            Item (row) number for the item in each of the itemParameter matrices to plot.
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric.
#' @param plotDensity       logical indicating if the focal distribution density should be plotted as a density curve (TRUE) or if it should be represented as an area gradient (FALSE). Defaults to gradient.
#' @param focalAbilities    If NULL, density is calculated theoretically from focal distribution. If not NULL, it must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution A string stating the distribution name to be used for density calculation. Only used if focalAbilities is NULL.
#' @param focalDistrExtra   Extra parameters for the focal group distribution function if needed.
#' @param from              value on the x-axis to serve as minimum for the plot
#' @param to                value on the x-axis to serve as maximum for the plot
#' @param thetaInt          value for the x-axis step for probabilities and density evaluation. Only used if focalAbilities is NULL.
#' @param colour            logical value indicating if the area gradient should be presented in colour when plotDensity is FALSE, or if the different lines should be presented in colour when plotDensity is TRUE.
#' @param highColour        character indicating the colour text name that should be used for high density regions.
#' @param main              text for plot main title.
#' @param xlab              text for x-axis label.
#' @param ylab              text for y-axis label.
#' @param iccText           text for legend title related to ICC curves.
#' @param focalIccText      legend for focal group ICC curve.
#' @param referenceIccText  legend for reference group ICC curve.
#' @param focalDensityText  legend for focal group density curve when plotDensity is TRUE. Text for legend title related to the colour gradient when plotDensity is FALSE.
#'
#' @return plotNCDIF  A ggplot object for the plot
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' # # Non Uniform - != guess DIF item close to focal distribution
#' PlotNcdif(iiItem = 46, itemParameters = dichotomousItemParameters, irtModel = "3pl",
#'           plotDensity = FALSE, main = "Item 46 Non uniform and different guessing DIF. 3PL")
#'
#' # # Non Uniform - != guess DIF item far from focal distribution
#' PlotNcdif(iiItem = 38, itemParameters = dichotomousItemParameters, irtModel = "3pl",
#'           plotDensity = FALSE, main = "Item 38 Non uniform and different guessing DIF. 3PL")
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
PlotNcdif <- function (iiItem, itemParameters, irtModel = "2pl", logistic = TRUE,

                       plotDensity = FALSE, focalAbilities = NULL, focalDistribution = "norm",
                       focalDistrExtra = list(mean = 0, sd = 1),
                       from = -5, to = 5, thetaInt = 0.01, colour = TRUE, highColour = "blue",
                       main = "", xlab = "Ability", ylab = "Probability",
                       iccText = "Group ICCs", focalIccText = "Focal group ICC",
                       referenceIccText = "Reference group ICC", focalDensityText = "Focal group density") {
#  require(ggplot2)

  # # Auxiliary functions
  CalculateProb <- switch(irtModel,
                          "1pl" = Calculate1plProb,
                          "2pl" = Calculate2plProb,
                          "3pl" = Calculate3plProb,
                          "grm" = CalculateGrmExp,
                          "pcm" = CalculatePcmExp,
                          stop("irtModel not known or not implemented"))

  FocalIcc <- function (x) {
    probability <- CalculateProb(thetaValue = x,
                                 itemParameters = matrix(itemParameters[["focal"]][iiItem, ], nrow = 1),
                                 logistic = logistic)
    return(probability)
  }

  ReferenceIcc <- function (x) {
    probability <- CalculateProb(thetaValue = x,
                                 itemParameters = matrix(itemParameters[["reference"]][iiItem, ], nrow = 1),
                                 logistic = logistic)
    return(probability)
  }

  # # Density calculation
  if (is.null(focalAbilities)) {
    thetaValue  <- seq(from = from, to = to, by = thetaInt)
    pars <- focalDistrExtra
    pars$x <-thetaValue
    focalWeight <- do.call(paste("d", focalDistribution, sep = ""), pars)
  } else {
    focalGroupDensity <- density(focalAbilities)
    thetaValue        <- focalGroupDensity$x
    focalWeight       <- focalGroupDensity$y
  }

  # # ICC information
  focalExpected     <- as.numeric(FocalIcc(thetaValue))
  referenceExpected <- as.numeric(ReferenceIcc(thetaValue))
  lowerExpected     <- apply(cbind(focalExpected, referenceExpected), 1, min)
  upperExpected     <- apply(cbind(focalExpected, referenceExpected), 1, max)

  plotData <- data.frame(thetaValue, lowerExpected, upperExpected, focalWeight)

  plotDataA   <- data.frame(thetaValue, icc = "Focal group ICC", expected = focalExpected)
  plotDataB   <- data.frame(thetaValue, icc = "Reference group ICC", expected = referenceExpected)
  plotDataICC <- rbind(plotDataA, plotDataB)

  # # Representation of density
  if (!plotDensity) {
    weightValues <- unique(focalWeight)

    plotNCDIF <- ggplot(plotData, aes(x = thetaValue))
    if (colour) {
      plotNCDIF <- plotNCDIF + geom_segment(aes(y = lowerExpected, xend = thetaValue, yend = upperExpected,
                                                colour = focalWeight))
    } else {
      plotNCDIF <- plotNCDIF + geom_segment(aes(y = lowerExpected, xend = thetaValue, yend = upperExpected,
                                                colour = focalWeight))
    }
  } else {
    plotNCDIF <- ggplot(plotDataICC, aes(x = thetaValue))
    plotDataC   <- data.frame(thetaValue, icc = "Focal group density", expected = focalWeight)
    plotDataICC <- rbind(plotDataICC, plotDataC)
  }

  # # Plotting of expected true scores by ability
  if (plotDensity & colour) {
    plotNCDIF <- plotNCDIF + geom_line(data = plotDataICC, aes(y = expected, colour = icc, linetype = icc), size = 0.8)
  } else {
    plotNCDIF <- plotNCDIF + geom_line(data = plotDataICC, aes(y = expected, linetype = icc), size = 0.8)
  }

  # # Legends
  if (plotDensity) {
    if (colour) {
      plotNCDIF <- plotNCDIF + scale_colour_hue(name  = iccText,
                                                breaks = c("Focal group ICC", "Reference group ICC", "Focal group density"),
                                                labels = c(focalIccText, referenceIccText, focalDensityText))
    }
    plotNCDIF <- plotNCDIF + scale_linetype(name   = iccText,
                                            breaks = c("Focal group ICC", "Reference group ICC", "Focal group density"),
                                            labels = c(focalIccText, referenceIccText, focalDensityText))
  } else {
    plotNCDIF <- plotNCDIF + scale_linetype(name   = iccText,
                                            breaks = c("Focal group ICC", "Reference group ICC"),
                                            labels = c(focalIccText, referenceIccText))
    if (colour) {
      highColour <- highColour
    } else {
      highColour = "black"
    }
    plotNCDIF <- plotNCDIF + scale_colour_gradient(name   = focalDensityText,
                                                   limits = c(0, max(focalWeight)),
                                                   low = "white",
                                                   high = highColour)
  }

  if (irtModel %in% c("grm", "pcm")) {
    if (ylab == "Probability") {
      ylab <- "Expected score"
    }
    ylimLow  <- 1
    ylimHigh <- ncol(itemParameters[["focal"]])
  } else {
    ylimLow  <- 0
    ylimHigh <- 1
  }

  plotNCDIF <- plotNCDIF + labs(title = main) + xlab(xlab) + ylab(ylab) + ylim(ylimLow, ylimHigh)

  return(plotNCDIF)
}
