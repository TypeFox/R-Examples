#' An example dataset
#'
#' Contains all the relevant information needed for the data reduction
#' some ARGUS-IV data from the University of Melbourne
#'
#' @name Melbourne
#' @docType data
#' @author David Philips \email{dphillip@@unimelb.edu.au}
#' @examples
#' data(Melbourne)
#' plotcorr(Melbourne$X)
NULL

#' The \code{redux} class
#'
#' An object class that is used throughout Ar-Ar_Redux
#' 
#' A list with the following items:
#' 
#' \code{labels}: a vector of strings denoting the names of the runs\cr
#' \code{num}: a vector of strings denoting the numerator isotopes\cr
#' \code{den}: a vector of strings denoting the denominator isotopes\cr
#' \code{intercepts}: a vector of logratio intercepts or values\cr
#' \code{covmat}: the covariance matrix of \code{intercepts}\cr
#' \code{irr}: a vector of strings denoting the irradiation runs\cr
#' \code{pos}: a vector of integers with the positions in the
#' irradiation stack \cr
#' \code{thedate}: a vector containing the acquisition dates and times\cr
#' \code{nlr}: a vector with the number of logratios per run\cr
#' \code{param}: a list of global parameters
#' @seealso param
#' @name redux
#' @docType class
NULL

#' The \code{timeresolved} class
#'
#' An object class containing time resolved multi-collector mass spectrometry data
#' 
#' A list with the following items:
#' 
#' \code{masses}: a vector of strings denoting the isotopes monitored in each run\cr
#' \code{irr}: a vector of strings denoting the irradiation runs\cr
#' \code{pos}: a vector of integers with the positions in the
#' irradiation stack \cr
#' \code{thedate}: a vector containing the acquisition dates and times\cr
#' \code{d}: a data table\cr
#' \code{thetime}: a matrix with the measurement times
#' @name timeresolved
#' @docType class
#' @seealso \code{\link{loaddata}}
NULL

#' The \code{PHdata} class
#'
#' An object class containing time resolved 'peak-hopping' mass
#' spectrometry data
#' 
#' A list with the following items:
#' 
#' \code{masses}: a vector of strings denoting the isotopes monitored in each run\cr
#' \code{signals}: a list with objects of class \code{\link{timeresolved}}, each corresponding
#' to a detector (i.e. length(\code{signals})==1 for single collector instruments).
#' @name PHdata
#' @docType class
#' @seealso \code{\link{loaddata}}
NULL

#' The \code{blankcorrected} class
#'
#' An object class containing blank-corrected mass spectrometry data
#'
#' Extends the class classes \code{\link{timeresolved}} and
#' \code{\link{PHdata}} by adding an additional list item
#' \code{blankindices} containg the index of the nearest
#' blank. \code{\link{fitlogratios}} uses this information to group
#' the samples during regression to 'time zero'.
#' @name blankcorrected
#' @docType class
NULL

#' The \code{logratios} class
#'
#' An object class containing logratio intercepts
#' 
#' A list with the following items:
#' 
#' \code{labels}: a vector of strings denoting the names of the runs\cr
#' \code{num}: a vector of strings denoting the numerator isotopes\cr
#' \code{den}: a vector of strings denoting the denominator isotopes\cr
#' \code{intercepts}: a vector of logratio intercepts or values\cr
#' \code{covmat}: the covariance matrix of \code{intercepts}\cr
#' \code{irr}: a vector of strings denoting the irradiation runs\cr
#' \code{pos}: a vector of integers with the positions in the
#' irradiation stack \cr
#' \code{thedate}: a vector containing the acquisition dates and times\cr
#' \code{nlr}: a vector with the number of logratios per run\cr
#' @name logratios
#' @docType class
NULL

#' The \code{results} class
#'
#' A list with the following items:
#' 
#' \code{labels}: a vector of strings denoting the names of the runs\cr
#' \code{intercepts}: a vector of ages\cr
#' \code{covmat}: the covariance matrix of \code{intercepts}\cr
#' \code{thedate}: a vector containing the acquisition dates and times\cr
#' @name results
#' @docType class
NULL
