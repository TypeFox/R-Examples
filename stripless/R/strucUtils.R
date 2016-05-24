#' @name strucUtils
#' @title Utility Functions for stripless
#' @description Unexported functions not intended to be directly called by users.
#'
#' @details \strong{Brief function descriptions:}
#' \describe{
#'  \item{check2lvl}{Checks for a 2 level design with center point}
#'  \item{make_2level_with_center}{Makes a 2 level design with center point to
#'   compare to a fraction of a \eqn{3^n} design}
#'   \item{makeSpacings}{Constructs a vector of spacings values for the "between"
#'   argument of \code{\link[lattice]{xyplot}}}
#'   }
#' @param spacings An integer vector
#' @param levelLength An integer vector giving the number of levels per factor
#'   at each level of the plotting hierarchy
#' @param d Condition list from strucParseFormula to be checked to see if it
#'   represents a 2 level design with center point.
#' @param sep sep argument for \code{paste()}.
#' @param lvls A list of length 3 character vectors that are the level names of
#'   factors in a design to be checked.
NULL
#
#' @rdname strucUtils
check2lvl <- function(d
  ,sep="." #sep character for paste()
  )
{
  lvl <- lapply(d,levels)
  all(sapply(lvl,length)==3) &&
  all(do.call(paste,c(d,list(sep=sep))) %in%
    make_2level_with_center(lvl,sep=sep))
}
#
#' @rdname strucUtils
make_2level_with_center <- function(
  lvls ## a list of length 3 character vectors that give factor levels)
  ,sep = "." ## sep character for paste
)
{
   corners <-do.call(paste,c(do.call(expand.grid,
   lapply(lvls,function(x)x[c(1,3)])),list(sep=sep)))
   center <- do.call(paste,c(lapply(lvls,function(x)x[2]),list(sep=sep)))
   half <- seq_len(2^(length(lvls)-1))
  c(corners[half],center,corners[-half])
}
#
#' @rdname strucUtils
makeSpacings <-   function(
  spacings
  ## vector giving spacing values. Will be extended as necessary by replicating last value
  ## Must all be >= 0
  ,levelLength
  ## vector giving sizes (lengths) of each level
)
{
  ## check vector values
  if(any(is.na(spacings))|| any(is.na(levelLength)))
    stop("NA values not allowed")
  stopifnot(all(spacings>=0),
            all(levelLength == floor(levelLength) & levelLength >=1))
  #####
  k <- length(levelLength)
  if(length(spacings)< k)
    spacings <- c(spacings,rep(spacings[length(spacings)],k))
  x <- NULL
  for(i in seq_len(k)) x <-
    c(x,rep(c(spacings[i],x),levelLength[i]-1))
  x
}

