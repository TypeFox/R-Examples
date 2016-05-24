#' Vectorized Version of outer
#' 
#' Vectorized \code{\link[base]{outer}}. 
#' 
#' @param x A \code{matrix}, \code{dataframe} or equal length \code{list} of 
#' vectors.
#' @param FUN A vectorized function.
#' @param \ldots Other arguments passed to the function supplied to \code{FUN}.
#' @return Returns a matrix with the vectorized \code{\link[base]{outer}} 
#' function.
#' @seealso \code{\link[base]{outer}},
#' \code{\link[stats]{cor}}
#' @author Vincent Zoonekynd, eddi of stackoverflow.com, and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @references \url{http://stackoverflow.com/a/9917425/1000343} \cr 
#' \url{http://stackoverflow.com/q/23817341/1000343}
#' @export
#' @rdname v_outer
#' @examples
#' #|------------------------------------------------------|
#' #|    SETTING UP VARIOUS FUNCTIONS THAT WILL BE USED    |
#' #|------------------------------------------------------|
#' pooled_sd <- function(x, y) {
#'     n1 <- length(x)
#'     n2 <- length(y)
#'     s1 <- sd(x)
#'     s2 <- sd(y)
#'     sqrt(((n1-1)*s1 + (n2-1)*s2)/((n1-1) + (n2-1)))
#' }
#' 
#' ## Effect Size: Cohen's d 
#' cohens_d <- function(x, y) {
#'     (mean(y) - mean(x))/pooled_sd(x, y)
#' }
#' 
#' 
#' ## Euclidean Distance
#' euc_dist <- function(x,y) sqrt(sum((x - y) ^ 2))
#' 
#' ## Cosine similarity
#' cos_sim <- function(x, y) x %*% y / sqrt(x%*%x * y%*%y)
#' 
#' sum2 <- function(x, y) sum(x, y)
#' arbitrary <- function(x, y) round(sqrt(sum(x)) - sum(y), digits=1)
#' #--------------------------------------------------------#
#' 
#' ## A data.frame
#' v_outer(mtcars, cor)
#' v_outer(mtcars, pooled_sd)
#' v_outer(mtcars[, 1:7], euc_dist)
#' v_outer(mtcars[, 1:7], sum2)
#' v_outer(mtcars[, 1:7], arbitrary)
#' 
#' ## mtcars as a list
#' mtcars2 <- lapply(mtcars[, 1:7], "[")
#' v_outer(mtcars2, cor)
#' v_outer(mtcars2, cor,  method = "spearman")
#' v_outer(mtcars2, pooled_sd)
#' v_outer(split(mtcars[["mpg"]], mtcars[["carb"]]), cohens_d)
#' v_outer(split(CO2[["uptake"]], CO2[["Plant"]]), cohens_d)
#' print(v_outer(mtcars[, 1:7], pooled_sd), digits = 1)
#' print(v_outer(mtcars[, 1:7], pooled_sd), digits = NULL)
#' v_outer(mtcars2, euc_dist)
#' v_outer(mtcars2, sum2)
#' v_outer(mtcars2, arbitrary)
#' 
#' ## A matrix
#' mat <- matrix(rbinom(500, 0:1, .45), ncol=10)
#' v_outer(mat, cos_sim)
#' v_outer(mat, euc_dist)
#' v_outer(mat, arbitrary)
#' 
#' \dontrun{
#' library(qdap)
#' wc3 <- function(x, y) sum(sapply(list(x, y), wc, byrow = FALSE))
#' L1 <- word_list(DATA$state, DATA$person)$cwl
#' (x <- v_outer(L1, wc3))
#' diag(x) <- (sapply(L1, length))
#' x
#' 
#' v_outer(with(DATA, wfm(state, person)), cos_sim)
#' with(DATA, Dissimilarity(state, person))
#' }
v_outer <-
function(x, FUN, ...){
	FUN
    UseMethod("v_outer")
}


#' @export
#' @method v_outer list
#' @rdname v_outer
v_outer.list <- 
function(x, FUN, ...){

    if (is.null(names(x))) {
        nms <- names(x) <- paste0("X", seq_along(x))
    } else {
        nms <- names(x)
    }

    z <- outer(
      nms, 
      nms, 
      Vectorize(function(i,j) FUN(unlist(x[[i]]), unlist(x[[j]]), ...))
    )
    dimnames(z) <- list(nms, nms)
    class(z) <- c("v_outer", class(z))
    z
}


#' @export
#' @method v_outer data.frame 
#' @rdname v_outer
v_outer.data.frame <- 
function(x, FUN, ...){

    nc <- ncol(x)
    mat <- matrix(rep(NA, nc^2), nc)
    for (i in 1:nc) {
        for (j in 1:nc) {
            mat[i, j] <- FUN(.subset2(x, i), .subset2(x, j))
        }
    }
    dimnames(mat) <- list(colnames(x), colnames(x))
    class(mat) <- c("v_outer", class(mat))
    mat
}


#' @export
#' @method v_outer matrix
#' @rdname v_outer
v_outer.matrix <- 
function(x, FUN, ...){

    x <- as.data.frame(x, stringsAsFactors = FALSE)

    nc <- ncol(x)
    mat <- matrix(rep(NA, nc^2), nc)
    for (i in 1:nc) {
        for (j in 1:nc) {
            mat[i, j] <- FUN(.subset2(x, i), .subset2(x, j))
        }
    }
    dimnames(mat) <- list(colnames(x), colnames(x))
    class(mat) <- c("v_outer", class(mat))
    mat
}

#' Prints a v_outer Object.
#' 
#' Prints a v_outer object.
#' 
#' @param x The v_outer object
#' @param digits Number of decimal places to print. 
#' @param \ldots ignored
#' @method print v_outer
#' @export
print.v_outer <-
function(x, digits = 3, ...) {
    WD <- options()[["width"]]
    options(width=3000)
    y <- unclass(x)
    if (is.numeric(y) & !is.null(digits)) {
        y <- round(y, digits = digits)
    }    
    print(y)
    options(width=WD)  
}
