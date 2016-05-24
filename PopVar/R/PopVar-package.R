#' Genomic breeding tools to 1) predict standard statistics and correlated response in plant populations, and 2) performs cross-validation to estimate genome-wide prediction accuracy
#'
#' @name PopVar-package
#' @aliases PopVar
#' @docType package
#' @description PopVar includes two functions useful for genome-based breeding: \itemize{
#'                \item \code{\link{pop.predict}} uses phenotypic and genotypic data from a set of individuals known as a training population (TP) and a set of candidate parents, which may or may not be included in the TP, to predict the mean (\eqn{\mu}), genetic variance (\emph{V_G}), and superior progeny value (\eqn{\mu}\emph{_sp}) of the half-diallel, or a defined set of pairwise bi-parental crosses between parents. When multiple traits are provided \code{pop.predict} will also predict the correlated responses and correlation between all pairwise traits. See \cite{Mohammadi, Tiede, and Smith (2015)} for further details.
#'                \item \code{\link{x.val}} performs cross-validation (CV) to estimate the accuracy of genome-wide prediction (otherwise known as genomic selection) for a specific training population (TP), i.e. a set of individuals for which phenotypic and genotypic data is available. Cross-validation can be conducted via one of two methods, see \code{Details} in \code{x.val} documentation for more information.
#'              } 
#'              The dataset \code{\link{think_barley.rda}}, previously described in \cite{Sallam et al. (2014)}, is provided as an example of the proper formatting of input files and also for users to become familiar with the functions within \code{PopVar}.
#' @author Tyler Tiede (maintainer) \email{tyler.tiede7@@gmail.com} and Mohsen Mohammadi
#'        
#'         Many thanks to Kevin Smith for supporting the project and Jeff Neyhart for helping with the initial 'debugging' efforts.
#' @import
#'    BGLR
#'    qtl
#'    rrBLUP
#' @references 
#'    Mohammadi M., T. Tiede, and K.P. Smith. 2015. PopVar: A genome-wide procedure for predicting genetic variance and correlated response in bi-parental breeding populations. Crop Sci. \emph{Accepted}.
#'    
#'    Sallam, A.H., J.B. Endelman, J-L. Jannink, and K.P. Smith. 2015. Assessing Genomic Selection Prediction Accuracy in a Dynamic Barley Breeding Population. Plant Gen. 8(1)

NULL