#' S3 Methods for Integration into \pkg{clue} Framework
#' 
#' These methods integrates the class \code{int.lsbclust} into the framework set out in
#' package \pkg{clue}.
#' 
#' @rdname lsbclusttoclue
#' @aliases cl_class_ids.int.lsbclust is.cl_partition.int.lsbclust is.cl_hard_partition.int.lsbclust
#' @param x An object of class \code{int.lsclust}
#' @method cl_class_ids int.lsbclust
#' @export cl_class_ids.int.lsbclust
cl_class_ids.int.lsbclust <- function(x) as.cl_class_ids(x$cluster)

#' @rdname lsbclusttoclue
#' @method is.cl_partition int.lsbclust
#' @export is.cl_partition.int.lsbclust
is.cl_partition.int.lsbclust <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @method is.cl_hard_partition int.lsbclust
#' @export is.cl_hard_partition.int.lsbclust
is.cl_hard_partition.int.lsbclust <- function(x) TRUE