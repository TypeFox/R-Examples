#' Small portion of RNA-seq data from plant physiology study.
#'
#'A subsample of the RNA-seq data from Baduel et al. studying Arabidopsis Arenosa physiology.
#'
#'@name baduel_small
#'@rdname baduel_small
#'@aliases baduel baduel_small baduel_gmt design expr_norm_corr
#'
#'@usage data(baduel_small)
#'
#'@references Baduel P, Arnold B, Weisman CM, Hunter B & Bomblies K, Habitat-Associated Life
#'History and Stress-Tolerance Variation in Arabidopsis Arenosa. \emph{Plant Physiology}: 01875, 2015.
#'
#'@references Agniel D, Hejblum BP, Variance component score test for
#' time-course gene set analysis of longitudinal RNA-seq data, \emph{submitted}, 2016.
#'
#'@format 3 objects\itemize{
#'\item{\code{design}:} a design matrix for the 48 measured samples, containing the following variables:\itemize{
#'  \item \code{SampleName} corresponding column names from \code{expr_norm_corr}
#'  \item \code{Population} a factor identifying the plant population
#'  \item \code{Age_weeks} numeric age of the plant at sampling time (in weeks)
#'  \item \code{Replicate} a purely technical variable as replicates are not from the same individual over weeks.
#'  Should not be used in analysis.
#'  \item \code{Verbnalized} a logical variable indicating whether the plant had undergone
#'  vernalization (exposition to cold and short day photoperiods)
#'}
#'\item{\code{baduel_gmt}:} a gmt file (see \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}})
#'\item{\code{expr_norm_corr}:} a numeric matrix
#'}
#'
#'@examples
#' \dontrun{
#' rm(list=ls())
#' data("baduel_small")
#' design$Intercept <- 1
#' design$PopulationKA <- 1*(design$Population=="KA")
#' design$AgeWeeks_Population <- design$Age_weeks*(design$Population=="KA")
#' design$Vernalized <- 1*design$Vernalized
#' design$Vernalized_Population <- design$Vernalized*(design$Population=="KA")
#'
#' set.seed(54321)
#' KAvsTBG <- tcgsa_seq(y=log2(expr_norm_corr+1), x=apply(as.matrix(design[, c("Intercept",
#'    "Vernalized", "Age_weeks", "Vernalized_Population", "AgeWeeks_Population"), drop=FALSE]),
#'        2, as.numeric),
#'                      phi=design[, c("PopulationKA"), drop=FALSE],
#'                      genesets=baduel_gmt$genesets[c(3,5)],
#'                      which_test = "permutation", which_weights = "loclin",
#'                      n_perm=5000, preprocessed = TRUE, doPlot = TRUE)
#'
#' set.seed(54321)
#' Cold <- tcgsa_seq(y=log2(expr_norm_corr+1), x=apply(as.matrix(design[, c("Intercept",
#'    "Age_weeks", "PopulationKA", "AgeWeeks_Population"), drop=FALSE]), 2, as.numeric),
#'                  phi=design[, c("Vernalized", "Vernalized_Population"), drop=FALSE],
#'                  genesets=baduel_gmt$genesets[c(3,5)],
#'                  which_test = "permutation", which_weights = "loclin",
#'                  n_perm=5000, preprocessed = TRUE, doPlot = TRUE)
#' }
#'
#'
#' @source \url{http://www.ncbi.nlm.nih.gov/bioproject/PRJNA312410}
#'
#' @keywords datasets
#' @docType data
NULL

