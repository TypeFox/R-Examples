##' Input checks and reference preparation for performance calculation
##'
##' Checks whether \code{r} and \code{p} are valid reference and predictions. If \code{p} is a
##' multiple of \code{r}, recycles \code{r} to the size and shape of \code{p}. If \code{r} has
##' additional length 1 dimensions (usually because dimensions were dropped from \code{p}), it is
##' shortend to the shape of \code{p}.
##'
##' In addition, any \code{NA}s in \code{p} are transferred to \code{r} so that these samples are
##' excluded from counting in \code{\link{nsamples}}.
##'
##' \code{checkrp} is automatically called by the performance functions, but doing so beforehand and
##' then setting \code{.checked = TRUE} can save time when several performance measures are to be
##' calculated on the same results.
##' @param r reference
##' @param p prediction
##' @return \code{r}, possibly recycled to length of \code{p} or with dimensions shortened to \code{p}.
##' @author Claudia Beleites
##' @export
##' @examples
##' ref <- softclassval:::ref
##' ref
##'
##' pred <- softclassval:::pred
##' pred
##' 
##' ref <- checkrp (r = ref, p = pred)
##' sens (r = ref, p = pred, .checked = TRUE)
checkrp <- function (r, p){
  if (is.null (dr <- dim (r))) dr <- length (r)
  if (is.null (dp <- dim (p))) dp <- length (p)

  recycle <- prod (dp) > prod (dr)
  
  if (prod (dr) > prod (dp))
    stop ("r must not be larger than p")

  dp <- dp [seq_along (dr)]
  
  if (dr [1] !=  dp [1])       
    stop ("numer of samples (nrow) of reference and prediction must be the same.")
  
  if (! is.na (dr [2]) && dr [2] != dp [2]) 
    stop ("p and r do not have the same number of columns (classes).")

  ## check class names if possible
  if (!is.null (colnames (r)) & ! is.null (colnames (p))){
  	reorder <- match (colnames (r), colnames (p))
  	
  	if (any (is.na (reorder)))
  		warning ("colnames of r (", paste (colnames (r)), ") and p  (", paste (colnames (p)), 
  					") do not match.")
  	else if (any (reorder != seq_len (dp [2])))
			warning ("Columns of r seem not to be in the same order as colnames of p.\n",
							 "To reorder them first, consider\n",
	             deparse (substitute (p)), " <- slice (", deparse (substitute (p)), 
							 ", j = c (", paste (reorder, collapse = ", "), "))")
  }
  	
    

  if (any (is.na (dp)) || any (dr != dp)) { # NA: p is shorter than r
   
    equaldims <- seq_len (min (which (is.na (dp) | dr != dp)) - 1) # first equal dims

    ## thereafter only length 1 dimensions are allowed
    if (any (dr [- equaldims] != 1)) 
      stop ("Dimension mismatch between r and p.")

    ## if p is shorter than r: shorten r
    if (any (is.na (dp))){
      a <- attributes (r)
      a$dim <- a$dim [equaldims]
      a$dimnames <- a$dimnames  [equaldims]
      mostattributes (r) <- a
    }
  }

  ## recycle if necessary
  if (recycle) {
    a <- attributes (r)
    r <- rep (r, length.out = length (p))
    mostattributes (r) <- attributes (p)
    dimnames (r) [seq_along (a$dimnames)] <- a$dimnames
  }

  ## make sure p == NA => r == NA
  is.na (r) <- is.na (p)
  
  r
}
.test (checkrp) <- function (){
  checkEquals (checkrp (ref,       pred                    ), ref      )
  checkEquals (checkrp (ref.array, pred.array              ), ref.array)
  checkEquals (checkrp (ref      , pred.array              ), ref.array, msg = "recycling r")
  checkEquals (checkrp (ref.array [,,1, drop = FALSE], pred), ref      , msg = "shortening r")

  checkException (checkrp (ref.array, pred                 )           , msg = "length (dim (r)) > length (dim (p))")
  checkException (checkrp (1 : 2,     1                    )           , msg = "nrow (r) != nrow (p)")
  checkException (checkrp (ref,       pred [, 1 : 2]       )           , msg = "ncol (r) != ncol (p)")
  
  tmp <- ref.array
  dim (tmp) <- c (dim(ref.array) [1 : 2], 1, dim (ref.array) [3])
  checkException (checkrp (tmp,       pred.array           )           , msg = "Dimension mismatch")

  ## check NAs are transferred correctly to reference
  tmp <- pred.array
  nas <- sample (length (pred.array), 10)
  tmp [nas] <- NA
  checkEquals (which (is.na (checkrp (ref,       tmp                    ))), sort (nas))
 
  ## warnings for colnames mismatches
  warnlevel <- options ()$warn
  options (warn = -1)
  checkEquals (checkrp (r = ref,      p = ref [, 3 : 1]       ), ref)
  options (warn = 2)
  checkException (checkrp (ref,       pred                    ))
  checkException (checkrp (r = ref,   p = ref [, 3 : 1]       ))
  options (warn = warnlevel)
}

##' Performance calculation for soft classification
##'
##' These performance measures can be used with prediction and reference being continuous class
##' memberships in [0, 1].
##'
##' The rows of \code{r} and \code{p} are considered the samples, columns will usually hold the
##' classes, and further dimensions are preserved but ignored.
##'
##' \code{r} must have the same number of rows and columns as \code{p}, all other dimensions may be
##' filled by recycling.
##'
##' \code{spec}, \code{ppv}, and \code{npv} use the symmetry between the performance measures as
##' described in the article and call \code{sens}.
##'
##' @rdname performance
##' @param r vector, matrix, or array with reference. 
##' @param p vector, matrix, or array with predictions
##' @param groups grouping variable for the averaging by \code{\link[base]{rowsum}}. If \code{NULL},
##' all samples (rows) are averaged.
##' @param operator the \code{\link[softclassval]{operators}} to be used 
##' @param drop should the results possibly be returned as vector instead of 1d array? (Note that
##' levels of \code{groups} are never dropped, you need to do that e.g. by
##' \code{\link[base]{factor}}.)
##' @param .checked for internal use: the inputs are guaranteed to be of same size and shape. If
##' \code{TRUE}, \code{confusion} omits input checking
##' @return numeric of size (ngroups x \code{dim (p) [-1]}) with the respective performance measure
##' @author Claudia Beleites
##' @seealso Operators: \code{\link{prd}}
##'
##' For the complete confusion matrix, \code{\link{confmat}}
##' @references see the literature in \code{citation ("softclassval")}
##' @export
##' @include softclassval.R
##' @examples
##'
##' ref <- softclassval:::ref
##' ref
##'
##' pred <- softclassval:::pred
##' pred
##'
##' ## Single elements or diagonal of confusion matrix
##' confusion (r = ref, p = pred)
confusion <- function (r = stop ("missing reference"), p = stop ("missing prediction"),
                       groups = NULL,
                       operator = "prd",
                       drop = FALSE, .checked = FALSE){
  operator <- match.fun (operator)
  if (! .checked)
    r <- checkrp (r, p)
  res <- operator (r = r, p = p)
  res <- groupsum (res, group = groups, dim = 1, reorder = FALSE, na.rm = TRUE)

  drop1d (res, drop = drop)
}
## testing by .test (sens)


##' Calculate the soft confusion matrix
##' 
##' @rdname performance
##' @export
##' @examples
##'
##' ## complete confusion matrix
##' cm <- confmat (r = softclassval:::ref, p = pred) [1,,]
##' cm
##' 
##' ## Sensitivity-Specificity matrix:
##' cm / rowSums (cm)
##'
##' ## Matrix with predictive values:
##' cm / rep (colSums (cm), each = nrow (cm))
confmat <- function (r = stop ("missing reference"), p = stop ("missing prediction"), ...){
	rx <- slice (r, j = rep (seq_len (ncol (r)), ncol (p)), drop = FALSE)
	colnames (rx) <- NULL
	
	px <- slice (p, j = rep (seq_len (ncol (p)), each = ncol (r)), drop = FALSE)
	colnames (px) <- NULL
	
	cm <- confusion (r = rx, p = px, ...)
	d <- dim (cm)
	dim (cm) <- c (d [1], ncol (r), ncol (p), d [- (1 : 2)])

	dn <- dimnames (p)
	if (is.null (dn)) dn <- rep (list (NULL), ndim (p))
	dn <- c (list (rownames (cm)), list (r = colnames (r)), dn [- 1L])
	names (dn) [3L] <- "p"
	dimnames (cm) <- dn
	cm
}
.test (confmat) <- function (){
  cm <- confmat (r = ref, p = pred)[1,,]
  warn <- options(warn = -1)$warn
  on.exit (options (warn = warn))
  for (r in colnames (ref))
    for (p in colnames (pred))
      checkEqualsNumeric (cm [r, p], confusion (r = ref [, r], p = pred [, p]))
  options (warn = warn)
  
  ## one sample only
  checkEquals (confmat (r = ref[1,,drop = FALSE], p = pred[1,,drop = FALSE])[1,,],
               structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0),
                         .Dim = c(3L, 3L),
                         .Dimnames = structure(list(r = c("A", "B", "C"), p = c("a", "b", "c")),
                           .Names = c("r", "p")))
               )
}


##' @param eps limit below which denominator is considered 0
##' @param op.dev does the operator measure deviation?
##' @param op.postproc if a post-processing function is needed after averaging, it can be given
##' here. See the example.
##' @rdname performance
##' @export
##' @examples
##'
##' ## sensitivities
##' sens (r = ref, p = pred)
sens <- function (r = stop ("missing reference"), p = stop ("missing prediction"), groups = NULL,
                  operator = "prd",
                  op.dev = dev (match.fun (operator)),
                  op.postproc = postproc (match.fun (operator)),
                  eps = 1e-8,
                  drop = FALSE,
                  .checked = FALSE){
  force (op.dev)
  force (op.postproc)

  if (! (isTRUE (op.dev) | isTRUE (! op.dev)))
    stop ("op.dev must either be TRUE or FALSE.")
  
  if (!is.null (op.postproc))
    POSTFUN <- match.fun (op.postproc)

  if (!.checked)
    r <- checkrp (r, p)                     # do the input checks.
  
  res <- confusion (r = r, p = p, groups = groups, operator = operator, drop = FALSE,
                    .checked = TRUE)
  nsmpl <- nsamples (r = r, groups = groups, operator = operator)

  if (any (nsmpl < res))
    warning ("denominator < enumerator.")
  nsmpl [nsmpl < eps] <- NA

  res <- res / nsmpl
  
  if (! is.null (op.postproc))          # e.g. root for wRMSE
    res <- POSTFUN (res)         
      
  if (op.dev)                           # for wMAE, wMSE, wRMSE, and the like
    res <- 1 - res

  res
}
.test (sens) <- function (){
  ops <- c ("strong", "weak", "prd", "and", "wMAE", "wMSE", "wRMSE")

  ## shape & names
  for (o in ops){
    ## vector
    tmp <- sens (r = v, p = v, operator = o)
    checkEquals (dim (tmp), 1L)
    checkTrue (is.null (dimnames (tmp))[[1]])
    checkTrue (is.null (names (tmp)))

    ## matrix
    tmp <- sens (r = v [1 : 4], p = m, operator = o)
    checkEquals (dim (tmp), c(1L, ncol (m)), msg = "matrix")
    checkEquals (dimnames (tmp), list (NULL, colnames (m)), msg = "matrix")
    checkTrue (is.null (names (tmp)), msg = "matrix")
    
    ## array
    tmp <- sens (r = rep (v [1 : 5], 2), p = pred.array, operator = o)
    checkEquals (dim (tmp), c (1, dim (pred.array) [-1]), msg = "array")
    checkEquals (dimnames (tmp), c (list (NULL), dimnames (pred.array) [-1]), msg = "array")
    checkTrue (is.null (names (tmp)), msg = "array")
  }
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array),
                      c (0.6, 0.32, NA, 0.85, 0.4, NA))
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "weak"),
                      c (0.675, 0.5, NA, 1, 1, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "strong"),
                      c (0.55, 0.2, NA, 0.75, 0, NA))
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "prd"),
                      c (0.6, 0.32, NA, 0.85, 0.4, NA))
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "and"),
                      c (0.2, NA, NA, 1, NA, NA))
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "wMAE"),
                      c (0.66, 0.68, NA, 1, 1, NA))
  
  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "wMSE"),
                      c (0.788, 0.86, NA, 1, 1, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, operator = "wRMSE"),
                      c (1 - sqrt (1.696/8), 1 - sqrt (.28/2), NA, 1, 1, NA))


  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups),
                      c (0.6, 0.6, NA, 0.32, NA, NA, 1, 0.6, NA, 0.4, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "weak"),
                      c (0.6, 0.8, NA, 0.5, NA, NA, 1, 1, NA, 1, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "strong"),
                      c (0.6, 1.4/3, NA, 0.2, NA, NA, 1, 1/3, NA, 0, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "prd"),
                      c (0.6, 0.6, NA, 0.32, NA, NA, 1, 0.6, NA, 0.4, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "and"),
                      c (0.2, NA, NA, NA, NA, NA, 1, NA, NA, NA, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "wMAE"),
                      c (0.6, 0.76, NA, 0.68, NA, NA, 1, 1, NA, 1, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "wMSE"),
                      c (0.728, 0.888, NA, 0.86, NA, NA, 1, 1, NA, 1, NA, NA))

  checkEqualsNumeric (sens (r = ref.array, p = pred.array, groups = ref.groups, operator = "wRMSE"),
                      1 - sqrt (1 - c (0.728, 0.888, NA, 0.86, NA, NA, 1, 1, NA, 1, NA, NA)))

}

##' @param ... handed to \code{sens}
##' @rdname performance
##' @export 
##' @examples
##'
##' ## specificities
##' spec (r = ref, p = pred)
spec <- function (r = stop ("missing reference"), p = stop ("missing prediction"), ...){
  sens (r = 1 - r, p = 1 - p, ...)
}

##' @rdname performance
##' @export
##' @examples
##'
##' ## predictive values
##' ppv (r = ref, p = pred)
ppv <- function (r = stop ("missing reference"), p = stop ("missing prediction"), ...,
                 .checked = FALSE){
  if (! .checked)
    r <- checkrp (r, p) 
  sens (r = p, p = r, ..., .checked = TRUE)
}
.test (ppv) <- function (){
  checkEquals (ppv (r = ref, p = pred.array),
               sens (r = pred.array, p = ref.array))
}

##' @rdname performance
##' @export
##' @examples
##' npv (r = ref, p = pred)
npv <- function (r = stop ("missing reference"), p = stop ("missing prediction"), ...,
                 .checked = FALSE){
  if (! .checked)
    r <- checkrp (r, p)
  sens (r = 1 - p, p = 1 - r, ..., .checked = TRUE)
}

.test (npv) <- function (){
  checkEquals (npv (r = ref, p = pred.array),
               sens (r = 1 - pred.array, p = 1 - ref.array))
}
