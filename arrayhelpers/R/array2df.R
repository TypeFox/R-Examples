##' array2df: Convert multidimensional array into matrix or data.frame
##' The "wide-format" array is converted into a "long-format" \code{matrix} or
##' \code{data.frame}.
##' 
##' If the resulting \code{data.frame} is too large to fit in memory, a
##' \code{matrix} might help.
##' 
##' The main benefit of this function is that it uses matrices as long as
##' possible. This can give large advantages in terms of memory consumption.
##' 
##' @title array2df 
##' @param x \code{array}
##' @param levels \code{list} with the levels for the dimensions of \code{x}.
##' 
##' If \code{levels[[i]]} is \code{NULL} no column is produced for this factor.
##' 
##' If \code{levels[[i]]} is \code{NA}, the result column is a numeric with
##'   range from \code{1} to \code{dim (x)[i]}
##'
##' If  \code{levels[[i]]} is \code{TRUE}, the levels are taken from the dimnames.
##' 
##' \code{names(levels)} yield the resulting column names.
##' @param matrix If \code{TRUE}, a numeric \code{matrix} rather than a
##'   \code{data.frame} is returned.
##' @param label.x Name for the column containing the \code{x} values.
##' @param na.rm should rows where the value of \code{x} is \code{NA} be removed?
##' @return A data.frame or matrix with \code{prod (dim (x))} rows and \code{length (dim (x)) + 1}
##' columns. 
##' @author Claudia Beleites
##' @export
##' @include arrayhelpers.R
##' @keywords array manip
##' @seealso \code{\link[utils]{stack}}
##' @include arrayhelpers.R
##' @examples
##' a <- arrayhelpers:::a
##' a
##' array2df (a)
##' array2df (a, matrix = TRUE)
##' 
##' array2df (a, levels = list(NULL, x = NA, c = NULL), label.x = "value")
##' 
##' array2df (a, levels = list(NULL, x = TRUE, c = c ("foo", "bar")), label.x = "value")
##' 
##' summary (array2df (a,
##'                    levels = list(NULL, x = NA, c = c ("foo", "bar")),
##'                    label.x = "value"))
##' 
##' summary (array2df (a,
##'                    levels = list(NULL, x = NA, c = c ("foo", "bar")),
##'                    label.x = "value",
##'                    matrix = TRUE))
##' 

array2df <- function (x, levels, matrix = FALSE, 
    label.x = deparse (substitute (x)), na.rm = FALSE) {

  dims <- dim (x)
  dn <- dimnames (x)
  if (is.null (dn))
   dn <- lapply (dims, function (x) NULL)

  ## prepare levels if not given
  if (missing (levels)){
    levels <- dn
    levels [sapply (levels, is.null)] <- NA
  }
  
  ## TRUE as abbreviation for "use dimnames"
  i <- sapply (dn, is.null)
  dn [i] <- lapply (dims [i], seq_len)
  i <- sapply (levels, isTRUE)
  levels [i] <- dn [i]
  
  ## get colnames from dimnames of x
  colnames <- names (dn)
  if (! is.null (names (levels))){
    i <- nzchar (names (levels))
    colnames [i] <- names (levels) [i]
  }

  ## prepare unfolding of array
  if (length (levels) != length (dims)) 
    stop ("Levels must have as many elements as x has dimensions.")

  cprod <- c (1, cumprod (dims))
  rprod <- c (rev (cumprod (rev (dims))), 1) [-1]

  ## the dimensions to unfold
  idim <- seq_along (dims) [!sapply (levels, is.null)]

  ## start with a matrix
  df <- matrix (x, nrow = length (x), ncol = length (idim) + 1)

  for (d in seq_along (idim))
    df [, d + 1] <- rep (seq_len (dims [idim [d]]),
                         each = cprod [idim [d]], times = rprod [idim [d]])

  ## delete rows where x is NA?
  if (na.rm) 
    df <- df [! is.na (df [,1]), ]
  
  if (!matrix) {
    df <- as.data.frame (df)

    for (d in seq_along (idim)) {

      ## turn to factor?
      if (!all (is.na (levels [[idim [d]]]))) {
        df [, d + 1] <- factor (df [, d + 1], labels = levels [[idim [d]]])
      }
      
    }
    
  }

  if (!is.null (colnames))
    colnames (df) <- c (label.x, colnames [idim])
  else
    colnames (df) <- c (label.x, paste ("d", idim, sep = ""))
  
  df
}

.test (array2df) <- function (){

  ## no dimnames
  x <- array (1:24, 4:1)
  checkEquals (array2df (x),
               structure (list (x = 1:24,
                                d1 = rep (1:4, 6),
                                d2 = rep (rep (1:3, each = 4), 2),
                                d3 = rep (1:2, each = 12),
                                d4 = rep (1, 24)),
                                .Names = c("x", "d1", "d2", "d3", "d4"),
                                row.names = c(NA, -24L),
                                class = "data.frame")
               )               

  checkEquals (array2df (x, matrix = TRUE), 
               structure (c (1:24,
                             rep (1:4, 6),
                             rep (rep (1:3, each = 4), 2),
                             rep (1:2, each = 12),
                             rep (1, 24)),
                          .Dim = c(24L, 5L),
                          .Dimnames = list(NULL, c("x", "d1", "d2", "d3", "d4")))
               )               
               
  checkEquals (array2df (a), structure(list(a = 1:24, rows = structure(c(1L, 2L, 3L, 4L, 1L, 
    2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 
    2L, 3L, 4L), .Label = c("a", "b", "c", "d"), class = "factor"), 
    columns = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 
    3L, 3L, 3L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L
    ), .Label = c("A", "B", "C"), class = "factor"), d3 = structure(c(1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("1", "2"), class = "factor")), .Names = c("a", 
    "rows", "columns", "d3"), row.names = c(NA, -24L), class = "data.frame"))

  checkTrue (is.matrix (array2df (a, matrix = TRUE)))

  checkEquals (array2df (a, list (TRUE, NULL, NA)),
               structure(list(a = 1:24, rows = structure(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L,
               3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), .Label = c("a", "b", "c",
               "d"), class = "factor"), d3 = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L,
               2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)), .Names = c("a", "rows", "d3"), row.names
               = c(NA, -24L), class = "data.frame")
               )  
}
