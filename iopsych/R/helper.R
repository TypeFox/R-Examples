# This folder contains helper functions


# These functions validate input -----------------------------------------------


#' Is a matrix a correlation matrix?
#' 
#' @param r_mat A correlation matrix.
#' @return Nothing. Will return an error if \code{r_mat} is not a correlation matrix
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.isCorMat
.isCorMat <- function(r_mat) {
    mat <- unname(r_mat)
    cond <- as.logical(c(0,0,0))
    cond[1] <- isSymmetric.matrix(mat) #Is symetric matrix
    cond[2] <- prod(diag(mat)) == 1 #All ones on diagonal
    cond[3] <- sum(as.integer(mat)) == dim(mat)[1] #All between (1 & -1)
    if(all(cond) == FALSE) {stop("Invalid Correlation Matrix")}
}

#' Are y_col and x_col appropriate indexs for r_mat?
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @return Nothing. Will result in an error if test fails.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.checkIndex
.checkIndex <- function(r_mat, y_col, x_col) {
    .isCorMat(r_mat)
    if(length(y_col) != 1) {stop("y_col must be of length one")}
    if(max(c(y_col, x_col)) > dim(r_mat)[1]) {stop("Invalid columns selected")}
}


# These functions reshape input ------------------------------------------------


#' Finds rxx and rxy for a correlation matrix
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @return Matrix rxx, and vector rxy.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.indexMat
.indexMat <- function(r_mat, y_col, x_col) {
    rxx <- r_mat[x_col, x_col]
    ryy <- r_mat[y_col, y_col]
    rxy <- r_mat[y_col, x_col]
    out <- list(rxx, ryy, rxy)
    names(out) <- c("rxx", "ryy", "rxy")
    return(out)
}

#' Unpacks key vector
#' 
#' @param key_vec A key vector.
#' @return Returns (1) A list of indices (2) A list of standardized weights.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.unpack
.unpack <- function(key_vec) {
    key_id <- which(key_vec != 0)
    if(length(key_id) == 0) {stop("Empty Key")}
    key_wt <- (key_vec[key_id])
    out <- list(key_id, key_wt)
    names(out) <- c("key_id","key_wt")
    return(out)
}

#' Unpacks key matrix
#' 
#' Works like .unpack but accepts matrices rather than vectors.
#'
#' @param key_mat A matrix of keys. Each row is a key.
#' @return Returns (1) A list of indices (2) A list of standardized weights.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.unpackMat
.unpackMat <- function(key_mat) {
    key_id <- apply(key_mat, 1, function(x) which(x != 0) )
    key_wt <- sapply(1:nrow(key_mat), function(x) key_mat[x,key_id[[x]]])
    out <- list(key_id, key_wt)
    names(out) <- c("key_id", "key_wt")
    return(out)
}

#' Appends a new variable into a correlation matrix.
#' 
#' @param r_mat A correlation matrix.
#' @param r_vec A vector of correlations to append to r_mat.
#' @param lab A column name for r_vec.
#' @return A larger correlation matrix.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' #data(dls2007)
#' #dat <- dls2007
#' #rxx <- dat[1:4, 2:5]
#' #corAdd(rxx, c(.1,.1,.1,.1), lab="V5")
#' @keywords internal
#' @rdname internal.corAdd
.corAdd <- function(r_mat, r_vec, lab=""){
    #Validate Input
    .isCorMat(r_mat)
    if(nrow(r_mat) != length(r_vec)) { stop("Input dimensions do not match.") }
    #Expand Correlation matrix
    new_mat <- cbind(r_mat, r_vec)
    out <- (rbind(new_mat, c(r_vec, 1)))
    colnames(out) <- c(colnames(r_mat), lab)
    return(out)
}

