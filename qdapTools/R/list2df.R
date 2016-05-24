#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{list2df} - Convert a named list of vectors to a dataframe.
#' 
#' @param list.object A named \code{\link[base]{list}} of vectors..
#' @param col1 Name for column 1 (the vector elements if converting a list or 
#' the rownames if converting a matrix).
#' @param col2 Name for column 2 (the names of the vectors).
#' @return \code{list2df} - Returns a dataframe with two columns.
#' @keywords collapse list
#' @export
#' @seealso \code{\link[qdapTools]{mtabulate}}
#' @rdname list2df
#' @examples
#' lst1 <- list(x=c("foo", "bar"), y=1:5)
#' list2df(lst1)
#' 
#' lst2 <- list(a=c("hello", "everybody"), b = mtcars[1:6, 1])
#' list2df(lst2, "col 1", "col 2")
#' 
#' matrix2df(mtcars)
#' matrix2df(cor(mtcars))
#' matrix2df(matrix(1:9, ncol=3))
#' 
#' vect2df(1:10)
#' vect2df(c(table(mtcars[, "gear"])))
#' 
#' list_df2df(list(mtcars, mtcars))
#' 
#' L1 <- list(a=1:10, b=1:6, c=5:8)
#' list_vect2df(L1)
#' 
#' L2 <- list(
#'   months=setNames(1:12, month.abb),
#'   numbers=1:6, 
#'   states=setNames(factor(state.name[1:4]), state.abb[1:4])
#' )
#' 
#' list_vect2df(L2)
#' 
#' set.seed(10)
#' cnts <- data.frame(month=month.name,
#'     matrix(sample(0:2, 36, TRUE), ncol=3))
#' 
#' counts2list(cnts[, -1], cnts[, 1])
#' df2matrix(cnts)
#' counts2list(df2matrix(cnts))
#' counts2list(t(df2matrix(cnts)))
#' 
#' mat <- matrix(1:9, ncol=3)
#' matrix2long(mat)
#' matrix2long(mtcars)
#' 
#' \dontrun{
#' library(qdap)
#' term <- c("the ", "she", " wh")
#' (out <- with(raj.act.1,  termco(dialogue, person, term)))
#' x <- counts(out)
#' 
#' counts2list(x[, -c(1:2)], x[, 1])
#' }
#' 
#' vect2list(LETTERS[1:10])
#' vect2list(LETTERS[1:10], numbered.names = TRUE)
#' x <- setNames(LETTERS[1:4], paste0("Element_", 1:4))
#' vect2list(x)
#' vect2list(x, FALSE)
#' vect2list(x, FALSE, TRUE)
list2df <- function(list.object, col1 = "X1", col2 = "X2") {

    ## Make sure the vectors have names; if not use numbers
    if (is.null(names(list.object))){
        names(list.object) <- seq_along(list.object)
    }

    dat <- data.frame(x = unlist(list.object, ,FALSE), 
        y = rep(names(list.object), sapply(list.object, length)), 
        stringsAsFactors = FALSE, check.names=FALSE, row.names=NULL)
    colnames(dat) <- c(col1, col2)
    dat
}

#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{matrix2df} - Convert a matrix to a dataframe and convert the rownames 
#' to the first column.
#' 
#' @param matrix.object A matrix or simple_triplet_matrix object.
#' @rdname list2df
#' @return \code{matrix2df} - Returns a dataframe.
#' @export
matrix2df <- function(matrix.object, col1 = "var1") {

        ## Convert simple_triplet_matrix to a matrix
    if("simple_triplet_matrix" %in% class(matrix.object)){
        matrix.object <- as.matrix(matrix.object)
    }
    
    if (is.null(rownames(matrix.object))) {
        rownames(matrix.object) <- 1:nrow(matrix.object)
    }
    dat <- data.frame(rownames(matrix.object), matrix.object, row.names=NULL, 
        stringsAsFactors = FALSE, check.names=FALSE)
    colnames(dat)[1] <- col1
    dat
}

#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{vect2df} - Convert a named vector to a dataframe.
#' 
#' @param vector.object A vector object.
#' @param order logical.  If \code{TRUE} the dataframe will be ordered.
#' @param rev logical. If \code{TRUE} and \code{order = TRUE} the dataframe will 
#' be ordered in descending order.
#' @rdname list2df
#' @return \code{vect2df} - Returns a dataframe.
#' @export
vect2df <- function(vector.object, col1 = "X1", col2 = "X2", order = TRUE, 
    rev = FALSE) {
    
    if (!is.vector(vector.object) | is.list(vector.object)){
        warning("Does not appear to be a vector: Results my be inconsistent")
    }
    if (is.null(names(vector.object))) {
        names(vector.object) <- paste0("x", pad(vector.object))
    }
    out <- data.frame(names(vector.object), vector.object, check.names=FALSE, 
    	stringsAsFactors = FALSE)
    colnames(out) <- c(col1, col2)
    if (order) {
        FUN <- match.fun(ifelse(rev, "rev", "c"))
        if (rev) {
            out <- out[order(-out[, col2]), ] 
        } else {
            out <- out[order(out[, col2]), ] 
        } 
        out[, col1] <- factor(out[, col1], levels=as.character(out[, col1]))
    }
    rownames(out) <- NULL
    out
}

#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{list_df2df} - Convert a list of equal numbered/named columns to a 
#' dataframe using the list names as the level two variable.
#' 
#' @param list.df.object A list of dataframes with equal number/named of columns.
#' @rdname list2df
#' @return \code{list_df2df} - Returns a dataframe.
#' @export
list_df2df <- function(list.df.object, col1 = "X1") {

    if (is.null(names(list.df.object))) {
        names(list.df.object) <- paste0("L", pad(1:length(list.df.object)))
    }
    list.names <- rep(names(list.df.object), sapply(list.df.object, nrow))
    out <- data.frame(list.names, do.call(rbind, list.df.object), 
        row.names=NULL, check.names=FALSE, stringsAsFactors = FALSE)
    colnames(out)[1] <- col1
    out
}

#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{list_vect2df} - Convert a list of named vectors to a hierarchical
#' dataframe.
#' 
#' @param list.vector.object A list of dataframes with equal number/named of 
#' columns.
#' @param col3 The name of the third column (\code{list_vect2df}).
#' @param \dots Further arguments passed to \code{vect2df}.
#' @rdname list2df
#' @return \code{list_vect2df} - Returns a dataframe.
#' @export
list_vect2df <- function(list.vector.object, col1 = "X1", col2 = "X2", 
    col3 = "X3", order=TRUE, ...) {

    fct <- sapply(list.vector.object, is.factor)
    if (sum(fct > 0)) {
        list.vector.object[fct] <- lapply(list.vector.object[fct], as.character)
    }
    cls <- class(unlist(list.vector.object))
    list.vector.object <- lapply(list.vector.object, methods::as, cls)

    list_df2df(lapply(list.vector.object, vect2df, col1=col2, col2=col3, 
    	order=order, ...), col1=col1)

}

#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{counts2list} - Convert a count matrix to a named list of elements.
#' 
#' @param mat A matrix of counts.
#' @param nm A character vector of names to assign to the list.
#' @rdname list2df
#' @return \code{counts2list} - Returns a list of elements.
#' @export
counts2list <- function(mat, nm = rownames(mat)) {
    nms <- colnames(mat)
    stats::setNames(apply(mat, 1, function(x) rep(nms, x)),  nm = nm)
}

#counts2list <- function(mat, nm = rownames(mat)) {
#    stats::setNames(lapply(1:nrow(mat), function(i) {
#        x <- unlist(mat[i, , drop = FALSE])
#        x <- x[x > 0]
#        rep(names(x), x)
#    }),  nm = nm)
#}


#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{vect2list} - Convert a vector to a named list.
#' 
#' @param use.names logical.  If \code{TRUE} and the vector is named, these 
#' names will be transferred to the list names.
#' @param numbered.names logical.  If \code{TRUE} padded numbers will be used
#' as list names.  If \code{FALSE} the vector elements themselves will become
#' the list names.
#' @rdname list2df
#' @return \code{vect2list} - Returns a list of named elements.
#' @export
vect2list <- function(vector.object, use.names = TRUE, numbered.names = FALSE){
    
    if (is.list(vector.object) | ! is.vector(vector.object)) {
        stop("`vector.object` is not a vector; results may be unstable")
    }

    if (!is.null(names(vector.object)) && use.names) {
        stats::setNames(as.list(vector.object), names(vector.object))
    } else {
        if (numbered.names) {
            stats::setNames(as.list(vector.object), pad(1:length(vector.object)))
        } else {
            stats::setNames(as.list(vector.object), as.character(vector.object))
        }
    }
}

#' List/Matrix/Vector to Dataframe/List/Matrix/Matrix
#' 
#' \code{df2matrix} - Convert a dataframe to a \code{matrix} and simultaneously 
#' move a column (default is the first column) to the rownames of a 
#' \code{matrix}.
#' 
#' @param data.frame.object A \code{data.frame} object.
#' @param i The column number or name to become the rownames of the 
#' \code{matrix}.
#' @rdname list2df
#' @return \code{df2matrix} - Returns a matrix.
#' @export
df2matrix <- function(data.frame.object, i = 1) {

    if (is.numeric(i)) {
        i <- colnames(data.frame.object)[i]
    }

    x <- as.matrix(data.frame.object[, !colnames(data.frame.object) %in% c(i)])
    row.names(x) <- data.frame.object[, i]
    x
}


#' List/Matrix/Vector to Dataframe/List/Matrix
#' 
#' \code{matrix2long} - Convert a matrix to a long format dataframe where column 
#' names become column 1, row names, column 2 and the values become column 3.
#' 
#' @rdname list2df
#' @return \code{matrix2long} - Returns a long format dataframe.
#' @export
matrix2long <- function(matrix.object, col1 = "cols", col2 = "rows", col3 = "vals"){

    if (is.null(rownames(matrix.object))) {
        rownames(matrix.object) <- seq_len(nrow(matrix.object))
    }

    if (is.null(colnames(matrix.object))) {
        colnames(matrix.object) <- seq_len(ncol(matrix.object))
    }

    out <- stats::setNames(data.frame(
        rep(colnames(matrix.object), each=nrow(matrix.object)),
        rep(rownames(matrix.object), ncol(matrix.object)),
        c(unlist(matrix.object)), 
        stringsAsFactors = FALSE
    ), c(col1, col2, col3))
    rownames(out) <- NULL
    out
}


