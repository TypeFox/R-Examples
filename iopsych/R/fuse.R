# These functions use fuse() to calculate correlations between item composites.


# Internal functions ----------------------------------------------------------


#' Correlations between item composites and the original correlation matrix
#' 
#' The key matrix is used to specify any number of weighted item composites.
#' The correlations between each specified composite and the original correlation
#' matrix are computed.
#'
#' @param r_mat A correlation matrix.
#' @param key_mat A matrix with one row for each composite and one column for 
#'        each item contained in r_mat. The value if each element corresponds
#'        to the weight given to an item.
#' @return A matrix of intercorrelations.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.fuseRmat
.fuseRmat <- function(r_mat, key_mat) {
    len <- dim(key_mat)[1]
    input <- .unpackMat(key_mat)
    fn <- function(x) {
        fuseVec(r_mat=r_mat, a=input$key_id[[x]], wt_a=input$key_wt[[x]])
    }
    return(sapply(seq(len), fn))
}
#' Computes the intercorrelations of item composites
#' 
#' The key matrix is used to specify any number of weighted item composites.
#' A correlation matrix of these composites is then computed and returned.
#'
#' @param r_mat A correlation matrix.
#' @param key_mat A matrix with one row for each composite and one column for 
#'        each item contained in r_mat. The value if each element corresponds
#'        to the weight given to an item.
#' @return A matrix of intercorrelations.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.fuseCom
.fuseCom <- function(r_mat, key_mat) {
    input <- .unpackMat(key_mat)
    unpackedfuse <- function(X, Y){
      ua <- .unpack(key_mat[X,])
      ub <- .unpack(key_mat[Y,])
      fuse(r_mat=r_mat, a=ua$key_id, b=ub$key_id, wt_a=ua$key_wt, wt_b=ub$key_wt)
    }
    len <- dim(key_mat)[1]
    return(outer(X=1:len, Y=1:len, FUN=Vectorize(unpackedfuse)))
}

#' The intercorrelation among items and composites made of these items.
#' 
#' The key matrix is used to specify any number of weighted item composites.
#' A correlation matrix of these composites and the original correlation matrix
#' is then computed and returned.
#'
#' @param r_mat A correlation matrix.
#' @param key_mat A matrix with one row for each composite and one column for 
#'        each item contained in r_mat. The value if each element corresponds
#'        to the weight given to an item.
#' @return A matrix of intercorrelations.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.fusefull
.fuseFull <- function(r_mat, key_mat) {
    swr <- .fuseRmat(r_mat=r_mat, key_mat=key_mat)
    sws <- .fuseCom(r_mat=r_mat, key_mat=key_mat)
    #Combine 4 matrices together.
    tp <- cbind(r_mat, swr)
    bt <- cbind(t(swr), sws)
    return(rbind(tp,bt))
}


# External functions ----------------------------------------------------------


#' Computes the correlation between two composites of items.
#' 
#' Computes the correlation between two composites of items. Composites may 
#' contain overalapping items. Items weights for each composite may be specified.
#'
#' @param r_mat A correlation matrix.
#' @param a The items used for composite A specified as a vector of column numbers.
#' @param b The items used for composite B specified as a vector of column numbers.
#' @param wt_a A vector containing the weights of each item in composite A.
#' @param wt_b A vector containing the weights of each item in composite B.
#' @return A correlation coefficient.
#' @author Allen Goebl and Jeff Jones
#' @references Lord, F.M. & Novick, M.R. (1968). \emph{Statisticl theories of 
#'             menal test scores.}, 97-98.
#' @examples
#' Rxx <- matrix(c(1.00, 0.25,  0.50,  0.61,
#'                 0.25, 1.00,  0.30,  0.10,
#'                 0.50, 0.30,  1.00, -0.30,
#'                 0.61, 0.10, -0.30,  1.00), 4, 4)
#' a   <- c(1, 3)
#' b   <- c(2, 4)
#' 
#' # Example using overlapping items and weights
#' Rxx  <- matrix(.3, 4, 4); diag(Rxx) <- 1
#' a    <- c(1, 2, 4)
#' b    <- c(2, 3)
#' wt_a <- c(.60, .25, .15)
#' wt_b <- c(2, 3)
#' 
#' fuse(r_mat = Rxx, a = a, b = b, wt_a = wt_a, wt_b = wt_b)
#' 
#' @export
fuse <- function(r_mat, a, b, wt_a=rep(1, length(a)), wt_b=rep(1, length(b))) {
    #Sanity check
    .isCorMat(r_mat)
    if(length(a) != length(wt_a)) {stop("wt_a must be the same length as a")}
    if(length(b) != length(wt_b)) {stop("wt_b must be the same length as b")}
    #Variances
    ab <- (wt_a %*% r_mat[a,b] %*% wt_b)
    aa <- (wt_a %*% r_mat[a,a] %*% wt_a)
    bb <- (wt_b %*% r_mat[b,b] %*% wt_b)
    #Equation
    return(ab / (sqrt(aa) * sqrt(bb)))
}

#' Computes the correlation between a composite and a vector of items.
#' 
#' @param r_mat A correlation matrix.
#' @param a The items used for composite A specified as a vector of column numbers.
#' @param wt_a A vector containing the weights of each item in composite A.
#' @param output Output can be set to "mat", to return a matrix made up of the
#' newly generated correlations appened to the original correlation matrix.
#' @return A vector of correlation coefficients.
#' @author Allen Goebl and Jeff Jones
#' @references Lord, F.M. & Novick, M.R. (1968). \emph{Statisticl theories of 
#'             mental test scores.}, 97-98.
#' @examples
#' data(dls2007)
#' dat <- dls2007
#' rxx <- dat[1:4, 2:5]
#' items <- c(1,3)
#' wt_a <- c(2,1)
#' 
#' fuseVec(r_mat=rxx, a=items)
#' fuseVec(r_mat=rxx, a=items, wt_a=wt_a, output="mat")
#' @export
fuseVec <- function(r_mat, a, wt_a=rep(1, length(a)), output="vec") {
    .isCorMat(r_mat)
    len <- nrow(r_mat)
    out.vec <- sapply(seq(len),
                  function(x) fuse(r_mat=r_mat, a=a, b=x, wt_a=wt_a, wt_b=1))
    if(output == "vec") { return(out.vec) }
    else if(output == "mat") { return(.corAdd(r_mat=r_mat, r_vec=out.vec)) }
    else (stop("The output type is incorrectly specified."))
}

#' The intercorrelation among items and composites made of these items.
#' 
#' The key matrix is used to specify any number of weighted item composites.
#' A correlation matrix of these composites and the original correlation matrix
#' is then computed and returned.
#'
#' @param r_mat A correlation matrix.
#' @param key_mat A matrix with one row for each composite and one column for 
#'        each item contained in r_mat. The value if each element corresponds
#'        to the weight given to an item.
#' @param type The type of output desired.
#' @return If type = cxc then a matrix of the intercorrelations between the specified 
#'         composites are returned. If type = cxr then the intercorrelations 
#'         between the original item and the specified composites are returned.
#'         If type = full then all the intercorrelations between both the original items and the specified composites are returned.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' Rxx <- matrix(c(1.00, 0.25,  0.50,  0.61,
#'                 0.25, 1.00,  0.30,  0.10,
#'                 0.50, 0.30,  1.00, -0.30,
#'                 0.61, 0.10, -0.30,  1.00), 4, 4); Rxx
#' 
#' # Single composite
#' Key <- matrix(c(1, 2, 3, -1), 1, 4); Key
#' 
#' fuseMat(r_mat = Rxx, key_mat = Key)
#' 
#' # Three composites
#' Key <- matrix(c(1, 2, 3, -1,
#'                 2, 1, 0, -2,
#'                 1, 1, 0,  0), 3, 4, byrow = TRUE)
#' 
#' fuseMat(Rxx, Key)
#' @export
fuseMat <- function(r_mat, key_mat, type="full"){
    #Check input
    .isCorMat(r_mat)
    if(nrow(r_mat) != ncol(key_mat)) { stop("key_mat does not match r_mat") }
    #Call specified function
    switch(type,
           full = { out <- .fuseFull(r_mat=r_mat, key_mat=key_mat)},
           rxc  = { out <- .fuseRmat(r_mat=r_mat, key_mat=key_mat)},
           cxc  = { out <- .fuseCom(r_mat=r_mat, key_mat=key_mat)})
    return(out)
}

