#' Differential Functioning of Items and Tests framework
#'
#' \code{DFIT} provides functions for calculating the differential item and test functioning
#' proposed by Raju et al. (1995).
#'
#' DFIT provides a set of functions to calculate the noncompensatory (NCDIF), compensatory (CDIF) and test level (DTF)
#' differential functioning indices for items and tests under Raju's (Raju, et al. 1995) DFIT framework.
#' It also provides functions for obtaining cut-off points for identifying differential functioning for these indices
#' following the Monte Carlo Item Parameter Replication approach proposed by Oshima et al. (2006).
#' 
#' This package also improves upon available DFIT software by allowing the covariance matrices for both focal and reference
#' groups to be used. This improves the obtained cut-off points, which result in type I error rates at the nominal level,
#' and increased power, when compared to the cut-off points obtained when using only the focal group item parameter
#' estimates and their estimate covariances (Cervantes, 2012). Furthermore, this package includes functions for obtaining
#' the asymptotic covariance matrices of item parameter estimates (currently only for dichotomous IRT models) and for
#' calculating the DFIT indices base on the focal group distribution as well as ability estimates for a sample from the
#' focal population are included; these enable ad hoc and a priori power calculations for given item parameters and sample
#' sizes to be possible with this package.
#' 
#' @references de Ayala, R. J., (2009). The theory and practice of item response theory. New York: The Guildford Press
#' @references Cervantes, V. H. (2012). On using the Item Parameter Replication (IPR) approach for power calculation of the noncompensatory differential item functioning (NCDIF) index (pp. 206-207). Proceedings of the V European Congress of Methodology. Santiago de Compostela, Spain: Universidade de Santiago de Compostela.
#' @references Cohen, A., Kim, S-H and Baker , F. (1993). Detection of differential item functioning in the Graded Response Moodel. Applied psychological measurement, 17(4), 335-350
#' @references Holland, P.W., and Thayer, D.T. (1988). Differential Item Performance and the Mantel-Haenszel Procedure. In H. Wainer and H.I. Braun (Eds.), Test Validity. Hillsdale, NJ: Erlbaum.
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically derived standard errors of Item Response Theory item parameter estimates. Journal of educational measurement, 41(2), 85--117.
#' @references Oshima, T. & Morris, S. (2008). Raju's Differential Functioning of Items and Tests (DFIT). Educational Measurement: Issues and Practice, 27(3), 43--50.
#' @references Oshima, T., Raju, N. & Nanda, A. (2006). A new method for assessing the statistical significance in the Differential Functioning of Items and Tests (DFIT) framework. Journal of educational measurement, 43(1), 1--17.
#' @references Raju, N. (1988). The area between two item characteristic cureves. Psychometricka, 53(4), 495--502.
#' @references Raju, N., Fortmann-Johnson, K., Kim, W., Morris, S., Nering, M. & Oshima, T. (2009). The item parameter replication method for detecting differential functioning in the polytomous DFIT framework. Applied psychological measurement, 33(2), 133--147.
#' @references Raju, N. S., van der Linden, W. J., & Fleer, P. F. (1995). An IRT-based internal measure of test bias with applications for differential item functioning. Applied Psychological Measurement, 19, 353--368.
#' @references Roussos, L., Schnipke, D. & Pashley, P. (1999). A generalized formula for the Mantel-Haenszel Differential Item Functioning parameter. Journal of educational and behavioral statistics, 24(3), 293--322.
#' @references Wright, K. (2011). Improvements for Differential Funtioning of Items and Tests (DFIT): Investigating the addition of reporting an effect size measure and power (Unpublished doctoral dissertation). Georgia State University, USA.
#'
#' @import simex mvtnorm ggplot2
#' @name DFIT
#' @docType package
NULL
