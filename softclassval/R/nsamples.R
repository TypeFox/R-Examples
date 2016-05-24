##' Count number of samples 
##'
##' Basically, the reference is summed up. For hard operators, the reference is hardened first: soft
##' values, i.e. \code{r} in (0, 1) are set to NA. 
##' @title Number of samples
##' @param r reference class labels with samples in rows.
##' @param groups grouping variable for the averaging by \code{\link[base]{rowsum}}. If \code{NULL},
##' all samples (rows) are averaged.
##' @param operator  the \code{\link[softclassval:operators]{operator}} to be used 
##' @param hard.operator optional: a logical determining whether only hard samples should be counted
##' @return number of samples in each group (rows) for each class (columns) and all further
##' dimensions of ref.
##' @author Claudia Beleites
##' @export
##' @include make01.R
##' @examples
##' ref <- softclassval:::ref
##' ref
##' nsamples (ref)
##' nsamples (ref, hard.operator = TRUE)
nsamples <- function  (r = r, groups = NULL, operator = "prd", hard.operator) {
  if (missing (hard.operator))
    hard.operator <- hard (match.fun (operator))
  else if (!is.logical (hard.operator) | is.na (hard.operator))
    stop ("hard.operator needs to be logical")

  if (hard.operator)
    r <- .make01 (r)
  
  groupsum (r, group = groups, dim = 1, reorder = FALSE, na.rm = TRUE)
}

.test (nsamples) <- function () {
  
  checkEqualsNumeric (nsamples (ref                             ), c (8, 2, 0))
  checkEqualsNumeric (nsamples (ref,        hard.operator = TRUE), c (5, 0, 0))
  checkEqualsNumeric (nsamples (pred.array                      ), c (6, 3.2, 0.8, 8, 2, 0))
  checkEqualsNumeric (nsamples (pred.array, hard.operator = TRUE), c (2, 2, 0, 5, 0, 0))

  checkEqualsNumeric (nsamples (ref,        groups = ref.groups                      ), c (5, 3, 0, 2, 0, 0))
  checkEqualsNumeric (nsamples (ref,        groups = ref.groups, hard.operator = TRUE), c (5, 0, 0, 0, 0, 0))
  checkEqualsNumeric (nsamples (pred.array, groups = ref.groups                      ),
                      c (3, 3, 1.6, 1.6, 0.4, 0.4, 5, 3, 0, 2, 0, 0))
  checkEqualsNumeric (nsamples (pred.array, groups = ref.groups, hard.operator = TRUE),
                      c (1, 1, 1, 1, 0, 0, 5, 0, 0, 0, 0, 0))

  checkEquals (nsamples (ref.array, groups = ref.groups),
               structure(c(5, 3, 0, 2, 0, 0, 5, 3, 0, 2, 0, 0), .Dim = c(2L, 3L, 2L),
                         .Dimnames = list (levels (ref.groups), colnames (ref.array), c("1", "2"))))
}
