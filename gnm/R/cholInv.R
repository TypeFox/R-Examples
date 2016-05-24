#  Copyright (C) 2006, 2010 David Firth and Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

cholInv <- function (mat,
                     eliminate = numeric(0),
                     onlyFirstCol = FALSE,
                     onlyNonElim = FALSE)
{
    .Deprecated(msg = paste("'cholInv' is deprecated as it is no longer used ",
                "by gnm."))
    m <- nrow(mat)
    n <- ncol(mat)
    if (length(eliminate) == 0) { ## the basic routine, no eliminated submatrix
        if (!is.matrix(mat))
            stop("mat is not a matrix")
        Rownames <- rownames(mat)
        Colnames <- colnames(mat)
        result <- chol2inv(chol(mat))
        if (!is.null(Rownames))
            colnames(result) <- Rownames
        if (!is.null(Colnames))
            rownames(result) <- Colnames
        if (onlyFirstCol)
            result <- result[, 1, drop = FALSE]
        return(result)
    }
##  Now allow for the possibility of an eliminated submatrix
    if (m != n)
        stop("mat must be a symmetric matrix")
    n <- nrow(mat)
    elim <- 1:n %in% eliminate
    diag.indices <- (n * (0:(n - 1)) + 1:n)
    T <- mat[diag.indices[eliminate]]
    if (any(T == 0))
      stop("an eliminated submatrix must have all diagonal entries non-zero.")
    W <- mat[!elim, !elim, drop = FALSE]
    U <- mat[elim, !elim, drop = FALSE]
    Ti <- 1/T
    k <- length(T)
    Ti.U <- Ti * U
    V.Ti <- t(Ti.U)
    Qmat <- W - crossprod(Ti.U, U)
    Qi <- cholInv(Qmat)
    result <- matrix(NA, if (onlyNonElim)
        n - k
    else n, if (onlyFirstCol)
        1
    else if (onlyNonElim)
        n - k
    else n)
    cols.notElim <- if (onlyFirstCol)
        1
    else if (onlyNonElim)
        1:(n - k)
    else !elim
    rows.notElim <- if (onlyNonElim)
        1:(n - k)
    else !elim
    if (onlyFirstCol)
        Qi <- Qi[, 1, drop = FALSE]
    result[rows.notElim, cols.notElim] <- Qi
    if (!onlyNonElim) {
        temp <- -crossprod(Qi, V.Ti)
        result[elim, cols.notElim] <- t(temp)
    }
    if (!onlyFirstCol && !onlyNonElim) {
        result[!elim, elim] <- temp
        temp <- crossprod(V.Ti, Qi) %*% V.Ti
        diag.indices <- k * (0:(k - 1)) + 1:k
        temp[diag.indices] <- Ti + temp[diag.indices]
        result[elim, elim] <- temp
    }
    theNames <- colnames(mat)
    rownames(result) <- if (onlyNonElim)
        theNames[!elim]
    else theNames
    colnames(result) <- if (onlyFirstCol)
        theNames[!elim][1]
    else if (onlyNonElim)
        theNames[!elim]
    else theNames
    result
}
