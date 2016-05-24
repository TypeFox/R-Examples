##' Retrieve data from the iDigBio specimen data repository.
##' 
##' @section About:
##' 
##' ridigbio provides an interface to the iDigBio data API described here: 
##' \url{https://www.idigbio.org/wiki/index.php/IDigBio_API}. With this package
##' you can retrieve specimen and media records from the iDigBio data 
##' repository. The iDigBio portal \url{https://portal.idigbio.org/} uses the 
##' same API so you should be able to retrieve the same information as shown in
##' the portal.
##' 
##' iDigBio contains nearly 30 million data records on musuem specimens held at
##' United States institutions. It also holds nearly 5 million images of these
##' specimens.
##'
##' @section Getting Started:
##' 
##' The main function is \code{\link{idig_search_records}} and reviewing its
##' documenation first with \code{?idig_search_records} is recommended.
##'
##' @section Limitations:
##' 
##' This package does not yet provide an interface to the mapping or the 
##' download APIs.
##'
##' @section Citing:
##' 
##' To cite the ridigbio package in your work, please use the following format:
##' 
##' Michonneau F, Collins M, Chamberlain SA (2016). ridigbio: An interface to iDigBio's search API that allows downloading specimen records. R package version 0.3.2. https://github.com/iDigBio/ridigbio
##'
##' @docType package
##' @name ridigbio
##' @title Retrieve data from the iDigBio specimen data repository.
##' @author Francois Michonneau \email{francois.michonneau@@gmail.com}
##' @author Matthew Collins \email{mcollins@@acis.ufl.edu}
NULL