#' Example Occurrence and Taxonomic Datasets of the Graptolithina from the Paleobiology Database
#'
#' Example datasets consisting of (a) occurrence data and (b) taxonomic data
#' downloaded from the Paleobiology Database API for the Graptolithina.

#' @name graptPBDB
#' @rdname graptPBDB
#' @aliases graptPBDB graptOccPBDB graptTaxaPBDB

#' @details
#' This example PBDB data is included here for testing functions involving occurrence data and taxonomy
#' in \code{paleotree}.

#' @format 
#' The example occurrence dataset (\code{graptOccPBDB}) is a data.frame consisting of 5900 occurrences (rows) and 35 variables (columns).
#' The example taxonomy dataset (\code{graptTaxaPBDB}) is a data.frame consisting of 364 formal taxa (rows) and 53 variables (columns).
#' Variables are coded in the 'pbdb' vocabulary of the PBDB API v1.2.

#' @seealso
#' \code{\link{taxonSortPBDBocc}}, \code{\link{occData2timeList}}, \code{\link{makePBDBtaxonTree}}, \code{\link{plotOccData}}

#' @source 
#' See examples for the full R code used to obtain the data from the API.
#' You can find the Paleobiology Database at http://paleobiodb.org
#' 
#' The occurrence data was entered by (in order of relative portion) P. Novack-Gottshall, M. Krause, M. Foote,
#' A. Hendy, T. Hanson, M. Sommers and others. This same data was authorized mainly by A. Miller,
#' W. Kiessling, M. Foote, A. Hendy, S. Holland, J. Sepkoski and others.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' \dontrun{
#' 
#' #original code used to obtain this dataset on March 21st, 2015
#' 		# using version 1.2 of the Paleobiology Database API
#'
#' # (sorry, URLs removed as they lead to the PBDB test server...)
#' 
#' save(graptOccPBDB,graptTaxaPBDB,file="graptPBDB.rdata")
#'
#' }
#'
#' # load archived example data
#' data(graptPBDB)
#'
#' # let's visualize who entered the majority of the occurrence data
#' pie(sort(table(graptOccPBDB$enterer)))
#' # and now who authorized it
#' pie(sort(table(graptOccPBDB$authorizer)))
#' # I apologize for using pie charts.
#' 
#' # Let's look at age resolution of these occurrences
#' hist(graptOccPBDB$early_age-graptOccPBDB$late_age,
#'		main="Age Resolution of Occurrences", xlab="Ma")
#' 
#' #distribution of taxa among taxonomic ranks
#' table(graptTaxaPBDB$taxon_rank)
#' barplot(table(graptTaxaPBDB$taxon_rank))
#'
NULL