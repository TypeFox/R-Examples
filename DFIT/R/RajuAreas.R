################################################################################
# # RajuAreas.R
# # R Versions: 2.14.1
# #
# # Author(s): Victor H. Cervantes
# #
# # Raju's DIF Area Measures
# # Description: Functions for calculating Raju's DIF Area Measures
# #
# # Inputs: NULL
# #
# # Outputs: functions
# #
# # File history:
# #   20120412: Creation
# #   20140425: Documentation adjusted to work with roxygen2. References
# #   added to documentation
################################################################################




################################################################################
# # Function UnsignedArea: Calculates Raju's Unsigned Area Measure index
# # for an item with given item parameters of focal and reference
# # groups.
################################################################################

#' Calculates Raju's Unsigned Area Measure index for an item with given item parameters of focal and reference groups.
#'
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return uam A numeric matrix with the Unsigned Area Measure values for all the item parameter in each set of itemParameterList
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' uam3pl <- UnsignedArea(itemParameters = dichotomousItemParameters, irtModel = "3pl",
#'                        subdivisions = 5000, logistic = TRUE)
#'
#' @references Cohen, A., Kim, S-H and Baker , F. (1993). Detection of differential item functioning in the Graded Response Moodel. Applied psychological measurement, 17(4), 335-350
#' @references Raju, N. (1988). The area between two item characteristic cureves. Psychometricka, 53(4), 495--502.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
UnsignedArea <- function (itemParameters, irtModel = "2pl", subdivisions = 5000, logistic = TRUE) {

  PointUnsignedDifference <- function (x, itemParameters, irtModel, logistic, ...) {
    itemDifference <- CalculateItemDifferences(thetaValue = x, itemParameters = itemParameters, irtModel = irtModel,
                                               logistic = logistic)
    out <- abs(itemDifference)
    return(out)
  }

  UnsignedAreaEqualDiscriminations <- function (focalDifficulty, referenceDifficulty, commonGuessing) {
    out <- (1 - commonGuessing) * abs(referenceDifficulty - focalDifficulty)
    return(out)
  }

  UnsignedAreaUnequalDiscriminations <- function (focalDifficulty, referenceDifficulty,
                                                  focalDiscrimination, referenceDiscrimination,
                                                  commonGuessing, logistic) {

    if (logistic) {
      kD <- 1
    } else {
      kD <- 1.702
    }

    discriminationFactor <- 2 * (referenceDiscrimination - focalDiscrimination) /
                            (kD * referenceDiscrimination * focalDiscrimination)
    difficultyDifference <- referenceDifficulty - focalDifficulty
    discriminationDifficultyFactor <- (kD * focalDiscrimination * referenceDiscrimination * difficultyDifference) /
                                      (referenceDiscrimination - focalDiscrimination)
    out <- (1 - commonGuessing) * abs(discriminationFactor * log(1 + exp(discriminationDifficultyFactor)) - difficultyDifference)
    return(out)
  }

  nItems <- nrow(itemParameters[["focal"]])
  uam <- rep(Inf, nItems)

  if (irtModel %in% c("1pl", "2pl", "3pl")) {
    if (irtModel == "1pl") {
      uam <- UnsignedAreaEqualDiscriminations(referenceDifficulty = itemParameters[["reference"]],
                                              focalDifficulty = itemParameters[["focal"]],
                                              commonGuessing = 0)
    }
    if (irtModel == "2pl") {
      for (ii in 1:nItems) {
        isEqualDiscriminations <- itemParameters[["reference"]][ii, 1] == itemParameters[["focal"]][ii, 1]
        if (isEqualDiscriminations) {
          uam[ii] <- UnsignedAreaEqualDiscriminations(focalDifficulty = itemParameters[["focal"]][ii, 2],
                                                      referenceDifficulty = itemParameters[["reference"]][ii, 2],
                                                      commonGuessing = 0)
        } else {
          uam[ii] <- UnsignedAreaUnequalDiscriminations(focalDifficulty = itemParameters[["focal"]][ii, 2],
                                                        referenceDifficulty = itemParameters[["reference"]][ii, 2],
                                                        focalDiscrimination = itemParameters[["focal"]][ii, 1],
                                                        referenceDiscrimination = itemParameters[["reference"]][ii, 1],
                                                        commonGuessing = 0, logistic = logistic)
        }
      }
    }
    if (irtModel == "3pl") {
      for (ii in 1:nItems) {
        isEqualDiscriminations <- itemParameters[["reference"]][ii, 1] == itemParameters[["focal"]][ii, 1]
        isEqualGuessing        <- itemParameters[["reference"]][ii, 3] == itemParameters[["focal"]][ii, 3]
        if (isEqualGuessing) {
          if (isEqualDiscriminations) {
            uam[ii] <- UnsignedAreaEqualDiscriminations(focalDifficulty = itemParameters[["focal"]][ii, 2],
                                                        referenceDifficulty = itemParameters[["reference"]][ii, 2],
                                                        commonGuessing = itemParameters[["focal"]][ii, 3])
          } else {
            uam[ii] <- UnsignedAreaUnequalDiscriminations(focalDifficulty = itemParameters[["focal"]][ii, 2],
                                                          referenceDifficulty = itemParameters[["reference"]][ii, 2],
                                                          focalDiscrimination = itemParameters[["focal"]][ii, 1],
                                                          referenceDiscrimination = itemParameters[["reference"]][ii, 1],
                                                          commonGuessing = itemParameters[["focal"]][ii, 3],
                                                          logistic = logistic)
          }
        }
      }
    }
  } else {
    for (ii in 1:nItems) {
      uam[ii] <- integrate(f = PointUnsignedDifference, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                           itemParameters = list(focal     = matrix(itemParameters[["focal"]][ii, ], nrow = 1),
                                                 reference = matrix(itemParameters[["reference"]][ii, ], nrow = 1)),
                           irtModel = irtModel, logistic = logistic)$value
    }
  }

  return(uam)
}








################################################################################
# # Function SignedArea: Calculates Raju's Signed Area Measure index
# # for an item with given item parameters of focal and reference
# # groups.
################################################################################

#' Calculates Raju's Signed Area Measure index for an item with given item parameters of focal and reference groups.
#'
#' @param itemParameters    A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return sam A numeric matrix with the Signed Area Measure values for all the item parameter in each set of itemParameterList
#'
#' @export
#'
#' @examples
#'
#' data(dichotomousItemParameters)
#' sam3pl <- SignedArea(itemParameters = dichotomousItemParameters, irtModel = "3pl",
#'                      subdivisions = 5000, logistic = TRUE)
#'
#' @references Cohen, A., Kim, S-H and Baker , F. (1993). Detection of differential item functioning in the Graded Response Moodel. Applied psychological measurement, 17(4), 335-350
#' @references Raju, N. (1988). The area between two item characteristic cureves. Psychometricka, 53(4), 495--502.
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
SignedArea <- function (itemParameters, irtModel = "2pl", subdivisions = 5000, logistic = TRUE) {

  PointSignedDifference <- function (x, itemParameters, irtModel, logistic, ...) {
    itemDifference <- CalculateItemDifferences(thetaValue = x, itemParameters = itemParameters, irtModel = irtModel,
                                               logistic = logistic)
    out <- itemDifference
    return(out)
  }

  SignedAreaDichotomous <- function (focalDifficulty, referenceDifficulty, commonGuessing) {
    out <- (1 - commonGuessing) * (referenceDifficulty - focalDifficulty)
    return(out)
  }

  nItems <- nrow(itemParameters[["focal"]])
  sam <- rep(Inf, nItems)

  if (irtModel %in% c("1pl", "2pl", "3pl")) {
    if (irtModel == "1pl") {
      sam <- SignedAreaDichotomous(referenceDifficulty = itemParameters[["reference"]],
                                              focalDifficulty = itemParameters[["focal"]],
                                              commonGuessing = 0)
    }
    if (irtModel == "2pl") {
          sam <- SignedAreaDichotomous(focalDifficulty = itemParameters[["focal"]][, 2],
                                                      referenceDifficulty = itemParameters[["reference"]][, 2],
                                                      commonGuessing = 0)
    }
    if (irtModel == "3pl") {
      for (ii in 1:nItems) {
        isEqualGuessing        <- itemParameters[["reference"]][ii, 3] == itemParameters[["focal"]][ii, 3]
        if (isEqualGuessing) {
            sam[ii] <- SignedAreaDichotomous(focalDifficulty = itemParameters[["focal"]][ii, 2],
                                                        referenceDifficulty = itemParameters[["reference"]][ii, 2],
                                                        commonGuessing = itemParameters[["focal"]][ii, 3])
        }
      }
    }
  } else {
    for (ii in 1:nItems) {
      sam[ii] <- integrate(f = PointSignedDifference, subdivisions = subdivisions, lower = -Inf, upper = Inf,
                           itemParameters = list(focal     = matrix(itemParameters[["focal"]][ii, ], nrow = 1),
                                                 reference = matrix(itemParameters[["reference"]][ii, ], nrow = 1)),
                           irtModel = irtModel, logistic = logistic)$value
    }
  }

  return(sam)
}
