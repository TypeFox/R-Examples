################################################################################
# # IPR.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Item Parameter Replication (Oshima, 2006) functions
# # Description: Functions for generating item parameters through Monte Carlo simulations and multivariate normality of
# #              estimates assumption
# #
# # Inputs: None
# #
# # Outputs: Functions
# #
# # File history:
# #   20120425: Creation
# #   20140421: Documentation adjusted to work with roxygen2
################################################################################




################################################################################
# # Function Ipr: Item parameter replication
################################################################################

#' Item parameter replication
#' @description Generates a sample of item parameters assuming multivariate normality of estimates
#'
#' @param itemParameters  A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items.
#' @param itemCovariances A list containing "focal" and "reference" lists of matrices of covariance for item estimates.
#' @param nReplicates     A numeric value indicating the number of replications to perform
#'
#' @return itemParameters A list with item parameters for focal and reference groups
#'
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
#'
#' @import simex mvtnorm
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
Ipr <- function (itemParameters, itemCovariances, nReplicates = 5000) {

#  require(simex)
#  require(mvtnorm)

  # # Data check
  if (nrow(itemParameters[["focal"]]) == nrow(itemParameters[["reference"]])) {
    nItem <- nrow(itemParameters[["focal"]])
  }  else {
    stop("There must be the same number of items for both 'focal' and 'reference' groups")
  }
  if (length(itemCovariances[["focal"]]) == length(itemCovariances[["reference"]])) {
    nCovs <- length(itemCovariances[["focal"]])
  }  else {
    stop("There must be the same number of item covariance matrices for both 'focal' and 'reference' groups")
  }

  if (length(nReplicates) > 1 | !is.numeric(nReplicates)) {
    stop("nReplicates must be a single numeric value")
  }
  if (nItem != nCovs) {
    stop("The number of item parameter vectors must be equal to the number of covariance matrices for each group")
  }

  # # IPR for each group
  itemPars <- as.numeric(t(itemParameters[["focal"]]))
  itemCovs <- simex::diag.block(itemCovariances[["focal"]])

  itemParameterFocalIpr <- mvtnorm::rmvnorm(n = nReplicates, mean = itemPars, sigma = itemCovs, method = "chol")
  itemParameterFocalIpr <- tapply(itemParameterFocalIpr , row(itemParameterFocalIpr),
                                  function (x) matrix(x, nrow = nItem, byrow = TRUE))

  itemPars <- as.numeric(t(itemParameters[["reference"]]))
  itemCovs <- simex::diag.block(itemCovariances[["reference"]])

  itemParameterReferenceIpr <- mvtnorm::rmvnorm(n = nReplicates, mean = itemPars, sigma = itemCovs, method = "chol")
  itemParameterReferenceIpr <- tapply(itemParameterReferenceIpr, row(itemParameterReferenceIpr),
                                      function (x) matrix(x, nrow = nItem, byrow = TRUE))

  # # Join the lists for each replication
  itemParameterList <- list()
  for (ii in seq(length(itemParameterFocalIpr))) {
    itemParameterList[[ii]]                <- list()
    itemParameterList[[ii]][["focal"]]     <- itemParameterFocalIpr[[ii]]
    itemParameterList[[ii]][["reference"]] <- itemParameterReferenceIpr[[ii]]
  }

  return(itemParameterList)
}








################################################################################
# # Function IprNcdif: NCDIF for Item parameter replication
################################################################################

#' NCDIF for Item parameter replication
#' @description Calculates the NCDIF index on a list of item parameters such as those produced by the Ipr function
#'
#' @param itemParameterList A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalAbilities    If NULL, NCDIF is calculated by numerical integration of focal distribution. If not NULL, must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution A string stating the distribution name to be used for integrating. Only used if focalAbilities is NULL.
#' @param focalDistrExtra   A list stating the extra parameters needed by the focal distribution function.
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities is NULL.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return ncdif A numeric matrix with the NCDIF values for all the item parameter in each set of itemParameterList
#'
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#' # # threePlNcdifIpr <- IprNcdif(itemParameterList = threePlIpr, irtModel = '3pl', logistic = TRUE)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
IprNcdif <- function (itemParameterList, irtModel = "2pl", focalAbilities = NULL, focalDistribution = "norm",
                      subdivisions = 5000, logistic = TRUE, focalDistrExtra = list(mean = 0, sd = 1)) {

  ncdif <- sapply(itemParameterList, function (x) Ncdif(x, irtModel = irtModel, focalAbilities = focalAbilities,
                                                        focalDistribution = focalDistribution,
                                                        focalDistrExtra = focalDistrExtra,
                                                        subdivisions = subdivisions, logistic = logistic))

  if (nrow(itemParameterList[[1]][['focal']]) == 1) {
    ncdif <- matrix(ncdif, nrow = 1)
  }

  return(ncdif)
}








################################################################################
# # Function IprUam: Unsigned Area Measure for Item parameter
# # replication
################################################################################

#' Unsigned Area Measure for Item parameter replication
#' @description Calculates Raju's Unsigned Area Measure index on a list of item parameters such as those produced by the Ipr function
#'
#' @param itemParameterList A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return uam A numeric matrix with the Unsigned Area Measure values for all the item parameter in each set of itemParameterList
#'
#' @references Raju, N. (1988). The area between two item characteristic cureves. Psychometricka, 53(4), 495--502.
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#' # # threePlUamIpr <- IprUam(itemParameterList = threePlIpr, irtModel = '3pl', logistic = TRUE)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
IprUam <- function (itemParameterList, irtModel = "2pl", subdivisions = 5000, logistic = TRUE) {

  uam <- sapply(itemParameterList, function (x) UnsignedArea(x, irtModel = irtModel, subdivisions = subdivisions,
                                                             logistic = logistic))

  if (nrow(itemParameterList[[1]][['focal']]) == 1) {
    uam <- matrix(uam, nrow = 1)
  }

  return(uam)
}








################################################################################
# # Function IprSam: Signed Area Measure for Item parameter replication
################################################################################

#' Signed Area Measure for Item parameter replication
#' @description Calculates Raju's Signed Area Measure index on a list of item parameters such as those produced by the Ipr function
#'
#' @param itemParameterList A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items.
#' @param irtModel          A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param subdivisions      A numeric value indicating the number of subdivisions for numerical integration.
#' @param logistic          A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return sam A numeric matrix with the Signed Area Measure values for all the item parameter in each set of itemParameterList
#'
#' @references Raju, N. (1988). The area between two item characteristic curves. Psychometricka, 53(4), 495--502.
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#' # # threePlSamIpr <- IprSam(itemParameterList = threePlIpr, irtModel = '3pl', logistic = TRUE)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
IprSam <- function (itemParameterList, irtModel = "2pl", subdivisions = 5000, logistic = TRUE) {

  sam <- sapply(itemParameterList, function (x) SignedArea(x, irtModel = irtModel, subdivisions = subdivisions,
                                                           logistic = logistic))

  if (nrow(itemParameterList[[1]][['focal']]) == 1) {
    sam <- matrix(sam, nrow = 1)
  }

  return(sam)
}








################################################################################
# # Function IprMh: Mantel Haenszel for Item parameter replication
################################################################################

#' Mantel Haenszel for Item parameter replication
#' @description Calculates the Mantel-Haenszel theoretical parameter under IRT assumptions on a list of item parameters such as those produced by the Ipr function
#'
#' @param itemParameterList     A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items
#' @param irtModel              A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".
#' @param focalDistribution     A string stating the distribution assumed for the focal group.
#' @param focalDistrExtra       A list stating the extra parameters needed by the focal distribution function.
#' @param referenceDistribution A string stating the distribution assumed for the reference group.
#' @param referenceDistrExtra   A list stating the extra parameters needed by the reference distribution function.
#' @param groupRatio            A positive value indicating how many members of the reference group are expected for each member of the focal group.
#' @param subdivisions          A numeric value indicating the number of subdivisions for numerical integration.
#' @param logistic              A logical value stating if the IRT model will use the logistic or the normal metric.
#'                         Defaults to using the logistic metric by fixing the D constant to 1.
#'                         If FALSE the constant is set to 1.702 so that the normal metric is used.
#'
#' @return mh A numeric matrix with the Mantel Haenszel values for all the item parameter in each set of itemParameterList
#'
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
#' @references Roussos, L., Schnipke, D. & Pashley, P. (1999). A generalized formula for the Mantel-Haenszel Differential Item Functioning parameter. Journal of educational and behavioral statistics, 24(3), 293--322.
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#' # # threePlMhIpr <- IprMh(itemParameterList = threePlIpr, irtModel = '3pl', logistic = TRUE)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
IprMh <- function (itemParameterList, irtModel = "2pl", focalDistribution = "norm",
                   focalDistrExtra = list(mean = 0, sd = 1), referenceDistribution = "norm",
                   referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1,
                   subdivisions = 5000, logistic = TRUE) {

  mh <- sapply(itemParameterList, function (x) IrtMh(itemParameters = x, irtModel = irtModel, subdivisions = subdivisions,
                                                     focalDistribution = focalDistribution,
                                                     focalDistrExtra = focalDistrExtra,
                                                     referenceDistribution = referenceDistribution,
                                                     referenceDistrExtra = referenceDistrExtra,
                                                     groupRatio = groupRatio,
                                                     logistic = logistic))

  if (nrow(itemParameterList[[1]][['focal']]) == 1) {
    mh <- matrix(mh, nrow = 1)
  }

  return(mh)
}








################################################################################
# # Function CutoffIpr: Cut-off points for Ipr generated estimates
################################################################################

#' Cut-off points for Ipr generated estimates
#' @description Calculates a given quantile cut-off point for each item on the IPR estimated items statistics. This
#' function  may produce the cut-off points for the NCDIF index, Signed and Unsigned Area Measures and the Mantel-Haenszel
#' statistic based on the Monte Carlo Item parameter replication approach. The quantiles may be calculated directly on the
#' output from the IprNcdif, IprSam, IprUam, and IprMh functions; the may be calculated by obtaining the corresponding
#' statistics for the item parameters simulated under the IPR approach; or by obatining both the simulated item parameters
#' and the statistics based on the item parameter values and their corresponging covariance matrices for the parameter
#' estimates. In the latter case, the user may choose to obtain the IPR simulated item parameters based only on the focal
#' group's covariance matrix as proposed by Oshima et al. (2006), or both focal and reference groups' matrices as proposed
#' by Cervantes (2012).
#'
#' @param iprStatistics         A numeric matrix with the statistics obtained for the simulated IPR item parameters or a list containing all the elements of the output of this function. If not NULL they will be used for calculating the cut-off points.
#' @param quantiles             A numeric vector with the quantiles to be calculated.
#' @param statistic             A character indicating which statistic will the cut-off point will be obtained for. If iprStatistics are provided, it is up to the user to correctly especify this string for it will only be informative; otherwise, it will be used to identify the statistic to be calculated. Should be one of "ncdif", "sam", "uam" or "mh".
#' @param itemParameterList     A list where each element is a list containing "focal" and "reference" item Parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items. Not used if iprStatistics are not NULL. If itemParameterList is not NULL, the statistic indicated with the argument "statistic" will be obtained for the set of itemParameterList, the corresponding arguments may be provided.
#' @param irtModel              A string stating the irtModel to be used. Should be one of "1pl", "2pl", "3pl", "grm" or "pcm". Not used if iprStatistics are not NULL.
#' @param focalAbilities        Only used if statistic is "ncdif". If NULL, NCDIF is calculated by numerical integration of focal distribution. If not NULL, must be a numerical vector containing the abilities for the individuals in the focal group.
#' @param focalDistribution     A string stating the distribution assumed for the focal group. Not used if iprStatistics are not NULL.
#' @param focalDistrExtra       A list stating the extra parameters needed by the focal distribution function. Not used if iprStatistics are not NULL.
#' @param referenceDistribution A string stating the distribution assumed for the reference group. Not used if iprStatistics are not NULL.
#' @param referenceDistrExtra   A list stating the extra parameters needed by the reference distribution function. Not used if iprStatistics are not NULL.
#' @param groupRatio            A positive value indicating how many members of the reference group are expected for each member of the focal group. Only used if iprStatistics are NULL and statistic is "mh".
#' @param focalSampleSize       A positive integer indicating the size of the focal group. Only used if itemCovariances is 'asymptotic'. Defaults to NULL.
#' @param referenceSampleSize   A positive integer indicating the size of the reference group. Only used if itemCovariances is 'asymptotic'. Defaults to NULL.
#' @param subdivisions          A numeric value indicating the number of subdivisions for numerical integration. Only used if focalAbilities and iprStatistics are NULL.
#' @param logistic              A logical value stating if the IRT model will use the logistic or the normal metric. Defaults to using the logistic metric by fixing the D constant to 1. If FALSE the constant is set to 1.702 so that the normal metric is used.
#' @param itemParameters        A list containing "focal" and "reference" item parameters. Item parameters are assumed to be on the same scale. Item parameters for each group should me a matrix with nrow equal to the number of items. Only used if both iprStatistics and itemParameterList are NULL. If used an itemParameterList from applying the IPR procedure will be simulated and the "statistic" will be calculated.
#' @param itemCovariances       Either a list containing "focal" and "reference" lists of matrices of covariance for item estimates or the string "asymptotic". Defaults to NULL. Only used if iprStatistics and itemParameterList are NULL, in all other cases the itemCovariances element of the returned list is equal to what is provided as value for these arguments.
#' @param nullGroup             If different from NULL and itemParameterList is NULL, a string equal to 'focal' or 'reference' to indicate which set of item parameters from itemParameters should be taken for the null hypothesis. If equal to NULL, itemParameterList will be generated using the given itemParameters for both groups.
#' @param nReplicates           A numeric value indicating the number of replications to perform. Only used if iprStatistics and itemParameterList are NULL.
#'
#' @return cutoff A list containing: 'itemParameters', NULL if not provided as argument, 'itemCovariances', NULL if not provided as argument, 'itemParameterList', NULL unless calculated from 'itemParameters' or provided as argument, 'iprStatistics' the matrix of 'statistics' provided as argument or calculated from 'itemParameterList', 'statistic' for which the IPR approach is used according to the provided argument, 'quantiles' the vector or matrix of calculated quantiles for each item
#'
#' @references Cervantes, V. H. (2012). On using the Item Parameter Replication (IPR) approach for power calculation of the noncompensatory differential item functioning (NCDIF) index (pp. 206-207). Proceedings of the V European Congress of Methodology. Santiago de Compostela, Spain: Universidade de Santiago de Compostela.
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
#' @references Raju, N. (1988). The area between two item characteristic curves. Psychometricka, 53(4), 495--502.
#' @references Roussos, L., Schnipke, D. & Pashley, P. (1999). A generalized formula for the Mantel-Haenszel Differential Item Functioning parameter. Journal of educational and behavioral statistics, 24(3), 293--322.
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
#' # # threePlAse[['reference']] <- AseIrt(itemParameters = dichotomousItemParameters[['focal']],
#' # #                                     logistic = TRUE, sampleSize = 500, irtModel = '3pl')
#' # # threePlIprCutoff <- CutoffIpr(itemParameters = dichotomousItemParameters,
#' # #                               itemCovariances = threePlAse, nullGroup = 'focal',
#' # #                               nReplicates = 1000, statistic = 'ncdif', irtModel = '3pl')
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
#'
CutoffIpr <- function (iprStatistics = NULL, quantiles, statistic = "ncdif",
                       itemParameterList = NULL, irtModel = "2pl", focalAbilities = NULL,
                       focalDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
                       referenceDistribution = "norm", referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1,
                       subdivisions = 5000, logistic = TRUE, itemParameters = NULL, itemCovariances = NULL,
                       nullGroup = NULL, focalSampleSize = NULL, referenceSampleSize = NULL, nReplicates = 5000) {

  if (is.null(iprStatistics)) {
    if (is.null(itemParameterList)) {
      if (is.null(itemParameters) || is.null(itemCovariances)) {
        stop("If iprStatistics and itemParameterList are not provided, itemParameters and itemCovaiances must be provided")
      }
      nullParameters <- itemParameters
      if (!is.null(nullGroup)) {
        if (nullGroup == 'focal') {
          nullParameters[['reference']] <- nullParameters[['focal']]
        } else if (nullGroup == 'reference') {
          nullParameters[['focal']] <- nullParameters[['reference']]
        } else {
          stop('nullGroup must be either NULL, "focal", or "reference"')
        } 
      } 
      if (is.character(itemCovariances) && itemCovariances == 'asymptotic') {
        itemCovariances <- list()
        itemCovariances[['focal']] <- AseIrt(itemParameters = nullParameters[['focal']],
                                             distribution = focalDistribution,
                                             distributionParameters = focalDistrExtra,
                                             sampleSize = focalSampleSize,
                                             irtModel = irtModel, logistic = logistic,
                                             subdivisions = subdivisions)
        itemCovariances[['reference']] <- AseIrt(itemParameters = nullParameters[['reference']],
                                                 distribution = referenceDistribution,
                                                 distributionParameters = referenceDistrExtra,
                                                 sampleSize = referenceSampleSize,
                                                 irtModel = irtModel, logistic = logistic,
                                                 subdivisions = subdivisions)
      } 
      itemParameterList <- Ipr(itemParameters = nullParameters, itemCovariances = itemCovariances,
                               nReplicates = nReplicates)

      if (irtModel == "3pl") {
        itemParameterList <- Bound3PlIpr(itemParameterList)
      }
    }

    if (statistic == 'ncdif') {
      iprStatistics <- IprNcdif(itemParameterList = itemParameterList, irtModel = irtModel,
                                focalAbilities = focalAbilities, focalDistribution = focalDistribution,
                                focalDistrExtra = focalDistrExtra,
                                subdivisions = subdivisions, logistic = logistic)
    }
    if (statistic == 'sam') {
      iprStatistics <- IprSam(itemParameterList = itemParameterList, irtModel = irtModel,
                              subdivisions = subdivisions, logistic = logistic)
    }
    if (statistic == 'uam') {
      iprStatistics <- IprUam(itemParameterList = itemParameterList, irtModel = irtModel,
                              subdivisions = subdivisions, logistic = logistic)
    }
    if (statistic == 'mh') {
      iprStatistics <- IprMh(itemParameterList = itemParameterList, irtModel = irtModel,
                             focalDistribution = focalDistribution, focalDistrExtra = focalDistrExtra,
                             referenceDistribution = referenceDistribution, referenceDistrExtra = referenceDistrExtra,
                             groupRatio = groupRatio,
                             subdivisions = subdivisions, logistic = logistic)
    }
  }

  if (is.list(iprStatistics)) {
    if (!all(names(iprStatistics) %in% c("itemParameters", "itemCovariances", "itemParameterList", "iprStatistics",
                                        "statistic", "quantiles"))) {
      stop("'iprStatistics' is a list, but does not contain all elements of 'CutoffIpr' output")
    }
    iprStatistics$iprStatistics <-  apply(iprStatistics$iprStatistics, 1, quantile, quantiles)
    cutoff <- iprStatistics
  } else {
    percentiles <- apply(iprStatistics, 1, quantile, quantiles)

    cutoff <- list(itemParameters = itemParameters, itemCovariances = itemCovariances,
                   itemParameterList = itemParameterList, iprStatistics = iprStatistics,
                   statistic = statistic, quantiles = percentiles)
  }

  return(cutoff)

}








################################################################################
# # Function Bound3PlIpr: Takes item parameters frp, Ipr and forces guessing to lie between 0 and 1
################################################################################

#' Takes item parameters frp, Ipr and forces guessing to lie between 0 and 1
#' @description Makes all simulated guessing values from a 3PL model that are outside the [0, 1] interval to be 0 or 1.
#'
#' @param itemParameterList     A list where each element is a list containing "focal" and "reference" item Parameters from a 3PL model. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items
#'
#' @return itemParameterList     A list where each element is a list containing "focal" and "reference" item Parameters where guessing parameters outside the [0, 1] interval are changed by 0 and 1.
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
#' # # threePlIpr <- Ipr(itemParameters = dichotomousItemParameters, itemCovariances = threePlAse,
#' # #                   nReplicates = 1000)
#' # # threePlIpr <- Bound3PlIpr(threePlIpr)
#'
#' @author Victor H. Cervantes <vcervantes at icfes.gov.co> <vhcervantesb at unal.edu.co>
Bound3PlIpr <- function (itemParameterList) {

  Change2Zero <- function (x) {
    isBelow <- which(x[, 3] < 0)
    x[isBelow, 3] <- 0
    return(x)
  }

  Change2One <- function (x) {
    isAbove <- which(x[, 3] > 1)
    x[isAbove, 3] <- 1
    return(x)
  }

  nCol <- sapply(itemParameterList, function (x) lapply(x, ncol))

  if (!all(nCol == 3)) {
    stop('Not all item parameter sets have three columns. They might not come from a 3PL model')
  }

  itemParameterList <- lapply(itemParameterList,
                              function (x) lapply(x, Change2Zero))

  itemParameterList <- lapply(itemParameterList,
                              function (x) lapply(x, Change2One))

  return(itemParameterList)
}
