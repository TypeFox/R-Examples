#' @title Rasch Model Parameters with pairwise
#' @name pairwise-package
#' @aliases pairwise
#' @docType package
#' @S3method logLik pers
#' @S3method plot grm 
#' @S3method plot pair
#' @S3method plot pairSE
#' @S3method plot pers
#' @S3method plot rfa
#' @S3method summary grm
#' @S3method summary pair
#' @S3method summary pairSE
#' @S3method summary pairwise.item.fit
#' @S3method summary pairwise.person.fit
#' @S3method summary pers
#' @S3method summary rfa
#' @importFrom graphics abline axis hist layout lines matlines matplot mtext par plot points rect segments symbols text title
#' @importFrom methods is
#' @importFrom stats aggregate complete.cases cov median na.omit pchisq sd var rnorm
#' @description Performs the explicit calculation -- not estimation! -- of the Rasch item parameters for dichotomous and polytomous response formats using a pairwise comparison approach (Choppin, 1968, 1985). On the basis of the item parameters, person parameters (WLE) are calculated according to Warm's weighted likelihood approach (Warm, 1989). Item- and person fit statistics and several functions for plotting are available.
#'
#' @details 
#' In case of dichotomous answer formats the item parameter calculation for the Rasch Model (Rasch, 1960), is based on the construction of a pairwise comparison matrix M\emph{nij} with entries f\emph{ij} representing the number of respondents who got item \emph{i} right and item \emph{j} wrong according to Choppin's (1968, 1985) conditional pairwise algorithm. 
#' 
#' For the calculation of the item thresholds and difficulty in case of polytomous answer formats, according to the Partial Credit Model (Masters, 1982), a generalization of the pairwise comparison algorithm is used. The construction of the pairwise comparison matrix is therefore extended to the comparison of answer frequencies for each category of each item. In this case, the pairwise comparison matrix M\emph{nicjc} with entries f\emph{icjc} represents the number of respondents who answered to item \emph{i} in category \emph{c} and to item \emph{j} in category \emph{c-1} widening Choppin's (1968, 1985) conditional  pairwise algorithm to polytomous item response formats. 
#' Within R this algorithm is simply realized by matrix multiplication.
#' 
#' In general, for both polytomous and dichotomous response formats, the benefit in applying this algorithm lies in it's capability to return stable item parameter 'estimates' even when using data with a relative high amount of missing values, as long as the items are still proper linked together.
#' 
#' The recent version of the package 'pairwise' computes item parameters for dichotomous and polytomous item responses  -- and a mixture of both -- according the partial credit model using the function \code{\link{pair}}.   
#' 
#' Based on the explicit calculated item parameters for a dataset, the person parameters may thereupon be estimated using any estimation approach. The function \code{\link{pers}} implemented in the package uses Warm's weighted likelihood approach (WLE) for estimation of the person parameters (Warm, 1989). When assessing person characteristics (abilities) using (rotated) booklet designs an 'incidence' matrix should be used, giving the information if the respective item was in the booklet (coded 1) given to the person or not (coded 0). Such a matrix can be constructed (out of a booklet allocation table) using the function \code{\link{make.incidenz}}.
#' 
#' Item- and person fit statistics, see functions \code{\link{pairwise.item.fit}} and \code{\link{pairwise.person.fit}} respectively, are calculated based on the squared and standardized residuals of observed and the expected person-item matrix. The implemented procedures for calculating the fit indices are based on the formulas given in Wright & Masters, (1982, p. 100), with further clarification given at \code{http://www.rasch.org/rmt/rmt34e.htm}. 
#' 
#' Further investigation of item fit can be done by using the function \code{\link{ptbis}} for point biserial correlations. For a graphical representation of the item fit, the function \code{\link{gif}} for plotting empirical and model derived category probability curves, or the function \code{\link{esc}} for plotting expected (and empirical) score curves, can be used.
#' 
#' The function \code{\link{iff}} plots or returns values of the item information function and the function \code{\link{tff}} plots or returns values of the test information function.
#' 
#' To detect multidimensionality within a set of Items a rasch residual factor analysis proposed by Wright (1996) and further discussed by Linacre (1998) can be performed using the function \code{\link{rfa}}. 
#' 
#' For a 'heuristic' model check the function \code{\link{grm}} makes the basic calculations for the graphical model check for dicho- or polytomous item response formats. The corresponding S3 plotting method is \code{\link{plot.grm}}. 
#'
#' @author Joerg-Henrik Heine <jhheine@@googlemail.com>
#' @references 
#' Choppin, B. (1968). Item Bank using Samplefree Calibration. \emph{Nature, 219}(5156), 870-872.
#' @references
#' Choppin, B. (1985). A fully conditional estimation procedure for Rasch model parameters. \emph{Evaluation in Education, 9}(1), 29-42.
#' @references
#' Heine, J. H. & Tarnai, Ch. (2015). Pairwise Rasch model item parameter recovery under sparse data conditions. \emph{Psychological Test and Assessment Modeling, 57}(1), 3–36. \url{http://www.psychologie-aktuell.com/index.php?id=inhaltlesen&tx_ttnews[tt_news]=3811&tx_ttnews[backPid]=204&cHash=e656c5b555}
#' @references
#' Heine, J. H. & Tarnai, Ch. (2011). Item-Parameter Bestimmung im Rasch-Modell bei unterschiedlichen Datenausfallmechanismen. \emph{Referat im 17. Workshop 'Angewandte Klassifikationsanalyse'} [Item parameter determination in the Rasch model for different missing data mechanisms. Talk at 17. workshop 'Applied classification analysis'], Landhaus Rothenberge, Muenster, Germany 09.-11.11.2011
#' @references
#' Heine, J. H., Tarnai, Ch. & Hartmann, F. G. (2011). Eine Methode zur Parameterbestimmung im Rasch-Modell bei fehlenden Werten. \emph{Vortrag auf der 10. Tagung der Fachgruppe Methoden & Evaluation der DGPs.} [A method for parameter estimation in the Rasch model for missing values. Paper presented at the 10th Meeting of the Section Methods & Evaluation of DGPs.] Bamberg, Germany, 21.09.2011 - 23.09. 2011.
#' @references Heine, J. H., & Tarnai, Ch. (2013). Die Pairwise-Methode zur Parameterschätzung im ordinalen Rasch-Modell. \emph{Vortrag auf der 11. Tagung der Fachgruppe Methoden & Evaluation der DGPs.} [The pairwise method for parameter estimation in the ordinal Rasch model. Paper presented at the 11th Meeting of the Section Methods & Evaluation of DGPs.] Klagenfurt, Austria, 19.09.2013 -  21.09. 2013.
#' @references Linacre, J. M. (1998). Detecting multidimensionality: which residual data-type works best? \emph{Journal of outcome measurement, 2}, 266–283.
#' @references Masters, G. N. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika, 47}(2), 149-174.
#' @references Rasch, G. (1960). \emph{Probabilistic models for some intelligence and attainment tests.} Copenhagen: Danmarks pædagogiske Institut.
#' @references Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory. \emph{Psychometrika, 54}(3), 427–450.
#' @references Wright, B. D., & Masters, G. N. (1982). \emph{Rating Scale Analysis.} Chicago: MESA Press.
#' @references Wright, B. D. (1996). Comparing Rasch measurement and factor analysis. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 3}(1), 3–24.

NULL
