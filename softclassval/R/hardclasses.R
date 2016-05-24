##' Convert to hard class labels
##'
##' \code{hardclasses} converts the soft class labels in \code{x} into a factor with hard class memberships and
##' \code{NA} for soft samples.
##'
##' \code{harden} hardens the soft 
##' @param x matrix or array holding the class memberships
##' @param classdim dimension that holds the classes, default columns
##' @param soft.name level for soft samples 
##' @param tol tolerance: samples with membership >= 1 - tol are considered to be hard samples of the
##' respective class.
##' @param drop see \code{\link[arrayhelpers]{drop1d}}
##' @return factor array of shape \code{dim (x) [-classdim]}
##' @seealso \code{\link[softclassval]{factor2matrix}} for the inverse 
##' @author Claudia Beleites
##' @export 
hardclasses <- function (x, classdim = 2L, soft.name = NA, tol = 1e-5, drop = TRUE){
  if (ndim (x) == 0) {                 # vector
    warning ("Using hardclasses (cbind (x, 1 - x)) instead.")
    x <- cbind (x, 1 - x)
    colnames (x) <- 1 : 0
  }
  
  classdim <- numericindex (x = dim (x), i = classdim, n = names (dimnames (x)))
  x <- aperm (x, c(seq_len (ndim (x))[-classdim], classdim))
  x <- makeNd (x, -2)
  olddims <- attr (x, "old")[[1]]
  
  if (any (abs(1 - rowSums (x)) > tol, na.rm = TRUE))
    warning ("Found samples with total membership != 1")
  
  if (is.null (classes <- colnames (x)))
      classes <- paste ("class", seq_len (ncol (x)), sep = "")
  
  x <- x >= 1 - tol                     # looses attributes!

  cl <- apply (x, 1, function (x) match (TRUE, x))
  
  if (! is.na (soft.name)){
    classes <- c (classes, soft.name)
    cl [is.na (cl)] <- length (classes)
  }

  cl <- structure (cl,
                   .Label    = classes, class = "factor",
                   .Dim      =      head (olddims$dim,      -1),
                   .Dimnames = lon (head (olddims$dimnames, -1)))
  drop1d (cl, drop = drop)
}

.test (hardclasses) <- function (){
  checkEquals (hardclasses (pred),
               factor (rep (letters [c (1, 2, NA, NA, NA)], 2), levels = letters [1 : 3]))

  checkEquals (hardclasses (pred, drop = FALSE), ensuredim (hardclasses (pred)))

  tmp <- pred
  dim (tmp) <- c (5, 2, 3)
  checkEquals (hardclasses (tmp, 3),
         structure (c (1L, 2L, NA, NA, NA, 1L, 2L, NA, NA, NA), .Dim = c(5L, 2L),
                    .Label = c("class1", "class2", "class3"), class = "factor"))
  
  ## vectors
  warn <- options(warn = 2)$warn
  on.exit (options (warn = warn))
  checkException (hardclasses (pred [,1]))
  options(warn = -1)
  checkEquals (hardclasses (pred [, 1]),
               factor (rep (c ("1", "0", NA, NA, NA), 2), levels = c ("1", "0")))
  options (warn = warn)

  ## NAs: missing predictions
  pred [2:3,] <- NA
  checkEquals (hardclasses (pred),
               factor (letters [c (1, NA, NA, NA, NA, 1, 2, NA, NA, NA)],
                       levels = letters [1 : 3]))
  pred [1,1] <- NA
  checkEquals (hardclasses (pred),
               factor (letters [c (NA, NA, NA, NA, NA, 1, 2, NA, NA, NA)],
                       levels = letters [1 : 3]))
  pred [6,2] <- NA
  checkEquals (hardclasses (pred),
               factor (letters [c (NA, NA, NA, NA, NA, 1, 2, NA, NA, NA)],
                       levels = letters [1 : 3]))
}

##' Mark operator as hard measure
##'
##' The operators may work only for hard classes (see \code{\link[softclassval:operators]{and}}). \code{hard (op)
##' == TRUE} marks hard operators.
##' 
##' @param op the operator (function)
##' @return logical indicating the type of operator. \code{NULL} if the attribute is missing.
##' @author Claudia Beleites
##' @seealso \code{\link{sens}} \code{\link[softclassval:operators]{and}}
##' @export 
##' @include softclassval.R
##' @include make01.R
##'
##' @examples
##'
##' hard (and)
##' myop <- function (r, p) p * (r == 1)
##' hard (myop) <- TRUE
##' 

hard <- function (op)
  attr (op, "hard")

##' @usage hard (op) <- value
##' @rdname hard
##' @param value logical indicating the operator type
##' @export "hard<-"
"hard<-" <- function (op, value){
  stopifnot (is.logical (value), !is.na (value))

  attr (op, "hard") <- value

  op
}

.test (hard) <- function (){
  myop <- function (){}
  checkTrue (is.null (hard (myop)))
  hard (myop) <- TRUE
  checkTrue (hard (myop))
  hard (myop) <- FALSE
  checkTrue (!hard (myop))
  checkException (hard (myop) <- NULL)
  checkException (hard (myop) <- NA)
}


##' @rdname hardclasses
##' @param closed logical indicating whether the system should be treated as closed-world (i.e. all
##' memberships add to 1)
##' @export harden
##' @include unittestdata.R
##' @examples
##' softclassval:::pred
##' harden (softclassval:::pred)
##' harden (softclassval:::pred, closed = FALSE)
##'
##' ## classical threshold at 0.5
##' harden (softclassval:::pred, tol = 0.5)
##'
##' ## grey zone: NA for memberships between 0.25 and 0.75
##' harden (softclassval:::pred, tol = 0.25)
##'
##' ## threshold at 0.7 = 0.5 + 0.2:
##' harden (softclassval:::pred - 0.2, tol = 0.5)
##' harden (softclassval:::pred - 0.2, tol = 0.5, closed = FALSE)

harden <- function (x, classdim = 2L, tol = 1e-6, closed = TRUE){
	x <- .make01 (x, tol = tol)

   if (closed && dim (x) [classdim] > 1L){
     ## in closed-world classification, one NA should "infect"
     ## all other classes of the same case

     nas <- colSums (aperm (x, c (classdim, seq_len (ndim (x)) [-classdim])))
     nas <- which (is.na (nas) | nas == 0, arr.ind = TRUE)
     nas <- as.matrix (nas)

     if (length (nas) > 0L)
       for (i in seq_len (dim  (x)[classdim])){
         tmp <- cbind (nas [,seq_len (classdim - 1)],
                       i,
                       nas [,seq_len (ncol (nas) - classdim + 1) + classdim - 1])
         x [tmp] <- NA
       }
   }
		
	x
}

.test (harden) <- function (){
  checkEquals (harden (as.matrix (v)),
               structure(c(0, NA, NA, 1, NA),
                         .Dim = c(5L, 1L),
                         .Dimnames = list(c("a", "b", "c", "d", "e"), NULL)))
  checkEquals (harden (m),
               structure(c(1, NA, NA, NA, 0, NA, NA, NA, 0, NA, NA, NA),
                         .Dim = c(4L, 3L),
                         .Dimnames = list(c("a", "b", "c", "d"), c("A", "B", "C")))
               )

  checkEquals (harden (pred.array),
               structure(c(1, 0, NA, NA, NA, 1, 0, NA, NA, NA,
                           0, 1, NA, NA, NA, 0, 1, NA, NA, NA,
                           0, 0, NA, NA, NA, 0, 0, NA, NA, NA,
                           1, 1, 1, 1, 1, NA, NA, NA, NA, NA,
                           0, 0, 0, 0, 0, NA, NA, NA, NA, NA,
                           0, 0, 0, 0, 0, NA, NA, NA, NA, NA),
                         .Dim = c(10L, 3L, 2L),
                         .Dimnames = list(NULL, c("a", "b", "c"), c("1", "2")))
               )
  checkEquals (harden (pred.array, closed = FALSE),
               structure(c(1, 0, NA, NA, NA, 1, 0, NA, NA, NA,
                           0, 1, NA, NA, NA, 0, 1, NA, NA, NA,
                           0, 0, 0,  NA, NA, 0, 0, 0,  NA, NA,
                           1, 1, 1, 1, 1, NA, NA, NA, NA, NA,
                           0, 0, 0, 0, 0, NA, NA, NA, NA, NA,
                           0, 0, 0, 0, 0,  0,  0,  0,  0,  0),
                         .Dim = c(10L, 3L, 2L),
                         .Dimnames = list(NULL, c("a", "b", "c"), c("1", "2")))
               )

  checkEquals (harden (pred.array, classdim = 3L),
               structure(c( 1,  0, NA, NA, NA, NA, NA, NA, NA, NA,
                           NA,  1, NA, NA, NA, NA, NA, NA, NA, NA,
                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            1,  1, NA, NA, NA, NA, NA, NA, NA, NA,
                           NA,  0, NA, NA, NA, NA, NA, NA, NA, NA,
                           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                         .Dim = c(10L, 3L, 2L),
                         .Dimnames = list(NULL, c("a", "b", "c"), c("1", "2")))
               )
}
