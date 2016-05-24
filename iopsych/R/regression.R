# This file contains all the functions for calculating R2 and beta weights


### Internal Functions --------------------------------------------------------


#' Correlation between weighted predictor composite and criterion.
#' 
#' @param wt_vec A vector of predictor weights.
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor, criterion correlations.
#' @return A correlation coefficent.
#' @note This is a simpler, faster version of the formula used for fuse().
#' @author Allen Goebl Jeff Jones
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' rxx <- dat[1:4, 1:4]
#' rxy <- dat[1:4, 5]
#'
#' #.solveWt(wt_vec=c(1,1,1,1), rxx=rxx, rxy=rxy)
#' #.solveWt(wt_vec=c(1,2,3,4), rxx=rxx, rxy=rxy)
#' @keywords internal
#' @rdname internal.solveWt
.solveWt <- function(wt_vec, rxx, rxy) {
    numer <- (t(wt_vec) %*% rxy)
    denom <- sqrt(t(wt_vec) %*% rxx %*% wt_vec)
    return(numer / denom)
}

#' Correlation between weighted predictor composite and criterion.
#' 
#' @param wt A vector of predictor weights, or a matrix of predictor weights with
#'        one column per predictor and one row per case.
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy_list A list of rxy vectors.
#' @return A matrix of correlation coefficent with one row per weight vector
#'         and one column per rxy vector.
#' @author Allen Goebl Jeff Jones
#' @note This function should be merged with the fuse functions and replace the
#'       other .solvewt functions.
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' rxx <- dat[1:4, 1:4]
#' rxy1 <- dat[1:4, 5]
#' rxy2 <- dat[1:4, 6]
#' rxy_list <- list(rxy1, rxy2)
#' 
#' wt1 <- c(1,1,1,1)
#' wt2 <- c(1,2,3,4)
#' wt_mat <- rbind(wt1, wt2)
#'
#' #.solveWtExp(wt=wt_mat, rxx=rxx, rxy_list=rxy_list)
#' @keywords internal
#' @rdname internal.solveWtExp
.solveWtExp <- function(wt, rxx, rxy_list) {
    solveRxyList <- function(wt_vec) {
        numer <- vapply(rxy_list, function(rxy) wt_vec %*% rxy, numeric(1))
        denom <- sqrt(t(wt_vec) %*% rxx %*% wt_vec)
        return(numer / denom)
    }
    #Correctly handle both vectors and matrices
    if(is.vector(wt)) { return(solveRxyList(wt)) }
    return(apply(wt, 1, FUN=solveRxyList))
}


#' Correlation between weighted predictor composite and criterion.
#' 
#' @param wt A vector of predictor weights or a list of vectors.
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor, criterion correlations.
#' @return A correlation coefficent.
#' @author Allen Goebl Jeff Jones
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' rxx <- dat[1:4, 1:4]
#' rxy <- dat[1:4, 5]
#' wt1 <- c(1,1,1,1)
#' wt2 <- c(1,2,3,4)
#' wt_list <- list(wt1, wt2)
#'
#' #.solveWtVec(wt=wt1, rxx=rxx, rxy=rxy)
#' #.solveWtVec(wt=wt2, rxx=rxx, rxy=rxy)
#' #.solveWtVec(wt=wt_list, rxx=rxx, rxy=rxy)
#' @keywords internal
#' @rdname internal.solveWtVec
.solveWtVec <- function(wt, rxx, rxy) {
    if(!is.list(wt)) { wt <- list(wt) }
    return(vapply(wt, function(x) .solveWt(wt_vec=x, rxx=rxx, rxy=rxy), numeric(1)))
}

#' Find regression weights
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor, criterion correlations.
#' @return A vector of regression weights.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.rmatBeta
.rmatBeta <- function(rxx, rxy) qr.solve(rxx) %*% (rxy)

#' Find regression weights and R2
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor, criterion correlations.
#' @return R2 and Regression weights
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.rmatReg
.rmatReg <- function(rxx, rxy) {
    beta <- qr.solve(rxx) %*% (rxy)
    R2 <- t(beta) %*% rxy
    out <- (list(R2, beta))
    names(out) <- c("R2", "beta")
    return(out)
}

#' Partially evaluated regression
#' 
#' Returns a function for calculating beta weights which has been partially
#' evalauted with respect to rxx.
#'
#' @param rxx A matrix of predictor intercorrelations.
#' @return Partially evaluated regression function.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.rmatBetaPE
.rmatBetaPE <- function(rxx){
    rxx_qr <- qr.solve(rxx)
    return(function(rxy) (rxx_qr %*% (rxy)))
}


# External functions ----------------------------------------------------------


#' Regression
#' 
#' @param r_mat A correlation matrix.
#' @param y_col The column representing the criterion variable.
#' @param x_col A vector of columns representing predictor variables.
#' @return Regression beta weights and R2.
#' @author Allen Goebl and Jeff Jones
#' @examples
#'Rs <- matrix(c(1.0, 0.2,  0.3, 0.4, -0.4,
#'               0.2, 1.0,  0.5, 0.1,  0.1,
#'               0.3, 0.5,  1.0, 0.2, -0.3,
#'               0.4, 0.1,  0.2, 1.0,  0.4,
#'              -0.4, 0.1, -0.3, 0.4,  1.0), 5, 5)
#'ys <- 5
#'xs <- 1:4
#'
#'rmatReg(Rs, ys, xs)
#' @export
rmatReg <- function(r_mat, y_col, x_col) {
    #Check Input
    .checkIndex(r_mat=r_mat, y_col=y_col, x_col=x_col)
    #Call Function
    term <- .indexMat(r_mat=r_mat, y_col=y_col, x_col=x_col)
    out <- .rmatReg(rxx=term$rxx, rxy=term$rxy)
    #Format Output
      #NOTE: Maybe use class here
    return(out)
}

#' Find r given arbitrary predictor weights
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param wt A vector of predictor weights or a list of multiple vectors.
#' @return The correlation between the weighted predictor composite and criterion.
#' @note This uses a simpler, faster version of the same formula used for fuse().
#' @author Allen Goebl and Jeff Jones
#' @examples
#' library(iopsych)
#' #Get Data
#' data(dls2007)
#' r_mat <- dls2007[1:6, 2:7]
#' 
#' #Get weights
#' unit_wt <- c(1,1,1,1)
#' other_wt <- c(1,2,1,.5)
#' wt_list <- list(unit_wt, other_wt)
#'
#' #Solve
#' solveWt(r_mat=r_mat, y_col=6, x_col=1:4, wt=unit_wt)
#' solveWt(r_mat=r_mat, y_col=6, x_col=1:4, wt=other_wt)
#' solveWt(r_mat=r_mat, y_col=6, x_col=1:4, wt=wt_list)
#' @export
solveWt <- function(r_mat, y_col, x_col, wt) {
    #Check Input
    .checkIndex(r_mat=r_mat, y_col=y_col, x_col=x_col)
    #Call Function
    term <- .indexMat(r_mat=r_mat, y_col=y_col, x_col=x_col)
    out <- .solveWtVec(rxx=term$rxx, rxy=term$rxy, wt=wt)
    #Format Output
      #NOTE: Maybe use class here
    return(out)
}

#' Find R2 given arbitrary predictor weights
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param wt A vector of predictor weights or a list of multiple vectors.
#' @return Regression R2.
#' @note This just calls solveWt() and squares the output.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' library(iopsych)
#' #Get Data
#' data(dls2007)
#' r_mat <- dls2007[1:6, 2:7]
#' 
#' #Get weights
#' unit_wt <- c(1,1,1,1)
#' other_wt <- c(1,2,1,.5)
#' wt_list <- list(unit_wt, other_wt)
#'
#' #Solve
#' solveWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=unit_wt)
#' solveWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=other_wt)
#' solveWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=wt_list)
#' @export
solveWtR2 <- function(r_mat, y_col, x_col, wt) {
    r <- solveWt(r_mat=r_mat, y_col=y_col, x_col=x_col, wt=wt)
    r2 <- (r^2)
    return(r2)
}

#' Partially evaluated regression
#' 
#' Returns a function for calculating beta weights and R2 which has been 
#' partially evalauted with respect to rxx.
#'
#' @param rxx A matrix of predictor intercorrelations.
#' @return Partially evaluated regression function.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' Rxx <- matrix(c(1.00, 0.25, 0.40,
#'                 0.25, 1.00, 0.30,
#'                 0.40, 0.30, 1.00), 3, 3)
#' 
#' rmatRegPE(Rxx)
#' @export
rmatRegPE <- function(rxx){
    #Check input
    .isCorMat(rxx)
    #Function
    rxx_qr <- qr.solve(rxx)
    fn <- function(rxy){
      beta <- rxx_qr %*% (rxy)
      R2 <- t(beta) %*% rxy
      out <- (list(beta, R2))
      names(out) <- c("beta", "R2")
      return(out)
    }
    return(fn)
}
