#' S3 Methods for Integration into \pkg{clue} Framework
#' 
#' These methods integrate the class \code{cds} into the framwork set out in
#' package \pkg{clue}. Use can therefore by made of \code{\link{cl_agreement}} to
#' calculate concordance measures between different solutions.
#' 
#' @rdname cdstoclue
#' @aliases cl_class_ids.cds is.cl_partition.cds is.cl_hard_partition.cds
#' @param x An object of class \code{cds}
#' @method cl_class_ids cds
#' @export cl_class_ids.cds
cl_class_ids.cds <- function(x) clue::as.cl_class_ids(x$grp)

#' @rdname cdstoclue
#' @method is.cl_partition cds  
#' @export is.cl_partition.cds
is.cl_partition.cds <- function(x) TRUE

#' @rdname cdstoclue
#' @method is.cl_hard_partition cds
#' @export is.cl_hard_partition.cds
is.cl_hard_partition.cds <- function(x) TRUE

#' @rdname cdstoclue
#' @method cl_class_ids cdsdata
#' @export cl_class_ids.cdsdata
cl_class_ids.cdsdata <- function(x) clue::as.cl_class_ids(x$grp.rs)

#' @rdname cdstoclue
#' @method is.cl_partition cdsdata
#' @export is.cl_partition.cdsdata
is.cl_partition.cdsdata <- function(x) TRUE

#' @rdname cdstoclue
#' @method is.cl_hard_partition cdsdata
#' @export is.cl_hard_partition.cdsdata
is.cl_hard_partition.cdsdata <- function(x) TRUE
