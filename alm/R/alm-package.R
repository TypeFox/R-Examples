#' R client for the open source article level metrics Lagotto application.
#'
#' An R interface to the RESTful API from the open source article level metrics
#' software Lagotto, created by the Public Library of Science (PLOS). A number of 
#' publishers are using Lagotto, so you can drop in a different base URL to the 
#' functions in this package to get to not only PLOS data, but
#' data for Crossref, and more.
#'
#' Authenication was required, but has now been removed going foward in Lagotto. 
#' You only need API keys for a few of the data providers that are running old 
#' versions of Lagotto, and that will change as they upgrade their Lagotto sofware.
#'
#' @name alm-package
#' @aliases alm
#' @docType package
#' @title R client for the open source article level metrics Lagotto application.
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @keywords package
NULL

#' Defunct functions in alm
#'
#' \itemize{
#'  \item \code{\link{alm_pubmedid}}: Function removed, you can get this info using alm_ids
#'  \item \code{\link{alm_pubmedcentid}}: Function removed, you can get this info using alm_ids
#'  \item \code{\link{almupdated}}: Function removed, you can get this info using alm_ids
#'  \item \code{\link{almdateupdated}}: Function removed, you can get this info using alm_ids
#'  \item \code{\link{alm_totals}}: Function removed, you can get this info using alm_ids
#' }
#'
#' @name alm-defunct
NULL
