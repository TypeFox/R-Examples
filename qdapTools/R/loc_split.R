#' Split Data Forms at Specified Locations
#' 
#' Split data forms at specified integer locations.
#' 
#' @param x A data form (\code{list}, \code{vector}, \code{data.frame}, 
#' \code{matrix}).
#' @param locs A vector of integer locations to split at.  If \code{locs} 
#' contains the index 1, it will be silently dropped.
#' @param names Optional vector of names to give to the list elements.
#' @param \ldots Ignored.
#' @return Returns of list of data forms broken at the \code{locs}.
#' @note Two dimensional object will retain dimension (i.e., \code{drop = FALSE} 
#' is used). 
#' @export
#' @seealso \code{\link[qdapTools]{run_split}},
#' \code{\link[qdapTools]{split_vector}}
#' \url{https://github.com/trinker/loc_split_example} for practical usage.
#' @examples
#' ## character
#' loc_split(LETTERS, c(4, 10, 16))
#' loc_split(LETTERS, c(4, 10, 16), c("dog", "cat", "chicken", "rabbit"))
#' 
#' ## numeric
#' loc_split(1:100, c(33, 66))
#' 
#' ## factor
#' (p_chng <- head(1 + cumsum(rle(as.character(CO2[["Plant"]]))[[1]]), -1))
#' loc_split(CO2[["Plant"]], p_chng)
#' 
#' ## list
#' loc_split(as.list(LETTERS), c(4, 10, 16))
#' 
#' ## data.frame
#' (vs_change <- head(1 + cumsum(rle(as.character(mtcars[["vs"]]))[[1]]), -1))
#' loc_split(mtcars, vs_change)
#' 
#' ## matrix
#' (mat <- matrix(1:50, nrow=10))
#' loc_split(mat, c(3, 6, 10))
loc_split <-
function(x, locs, names = NULL, ...) {

    locs
    names
    UseMethod("loc_split")

}


#' @export
#' @method loc_split list
#' @rdname loc_split
loc_split.list <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_vector(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split data.frame
#' @rdname loc_split
loc_split.data.frame <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_mat(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split matrix 
#' @rdname loc_split
loc_split.matrix <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_mat(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split numeric 
#' @rdname loc_split
loc_split.numeric <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_vector(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split factor
#' @rdname loc_split
loc_split.factor <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_vector(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split character 
#' @rdname loc_split
loc_split.character <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_vector(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

#' @export
#' @method loc_split default 
#' @rdname loc_split
loc_split.default <-
function(x, locs, names = NULL, ...) {

    names <- name_len_check(locs, names)
    out <- loc_split_vector(x, locs, ...)
    if(!is.null(names)) names(out) <- names
    out
}

loc_split_vector <- function(x, locs){
    starts <- c(1, locs)
    Map(function(s, e) {x[s:e]}, starts, c(locs - 1, length(x)))
}

## loc_split_vector <- #old version retired 9/9/15
## function(x, locs, names = NULL, ...) {
## 
##     locs <- locs[!locs %in% "1"]		
##     if (length(x) < max(locs)) stop("One or more `locs` elements exceeds length of `x`")
##     stats::setNames(split(x, cut(seq_along(x), c(0, locs - 1, length(x)))) , NULL)
## 
## }

loc_split_mat <- 
function(x, locs, names = NULL, ...) {

    locs <- locs[!locs %in% "1"]		
    len <- nrow(x)
    if (len < max(locs)) stop("One or more `locs` elements exceeds nrow of `x`")

    starts <- c(1, locs)
    Map(function(s, e) {x[s:e, ,drop=FALSE]}, starts, c(locs - 1, nrow(x)))

}

## loc_split_mat <- #old version retired 9/9/15
## function(x, locs, names = NULL, ...) {
## 
##     locs <- locs[!locs %in% "1"]		
##     len <- nrow(x)
##     if (len < max(locs)) stop("One or more `locs` elements exceeds nrow of `x`")
##     seqs <- seq_len(len)
## 
##     splitseqs <- split(seqs, cut(seqs, c(0, locs - 1, len)))
##     stats::setNames(lapply(splitseqs, function(i) x[i, ,drop=FALSE]), NULL)
## 
## }


name_len_check <- function(locs, names) {
    
    if (is.null(names)) return(names)
    check <- length(locs) + 1 == length(names)
    if(!check) warning("length of `names` muse be equal to length of `locs` + 1; ignoring `names`")
    if (!check) NULL else names
}



