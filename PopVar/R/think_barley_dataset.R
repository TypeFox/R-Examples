#' An example barley dataset
#' @name think_barley.rda
#' @aliases G.in_ex
#'          y.in_ex
#'          map.in_ex
#'          cross.tab_ex
#' @docType data
#' @description A sample dataset, previously described in \cite{Sallam et al. (2014)} is provided as an example of the proper formatting of input files and also for users to become familiar with \code{PopVar}; the \code{think_barley} dataset is useful in demonstrating both \code{\link{pop.predict}} and \code{\link{x.val}}. Note that a number of entries are missing data for one or both traits,
#'         which is representative of a real breeding scenario where phenotypic data may not be available for all parent candidates.
#' @format  The names of the example files are: 
#'        \describe{
#'          \item{G.in_ex}{A set of 245 barley lines genotyped with 742 SNP markers}
#'          \item{y.in_ex}{Phenotypes of four traits for a portion of the 245 barley lines, Fusarium head blight (FHB), deoxynivalenol (DON) in ppm, grain yield in bushels/acre, and plant height in cm.}
#'          \item{map.in_ex}{Genetic map (i.e. chromosome assignment and genetic distance (cM) between markers) of the 742 SNP markers based on \cite{Munoz-Amatriain et al., 2011}}
#'          \item{cross.tab_ex}{A table of user-defined crosses}
#'          }
#' @references
#'  Sallam, A.H., J.B. Endelman, J-L. Jannink, and K.P. Smith. 2015. Assessing Genomic Selection Prediction Accuracy in a Dynamic Barley Breeding Population. Plant Gen. 8(1)

NULL