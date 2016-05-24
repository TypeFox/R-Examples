MRM <- function(formula = formula(data), data = sys.parent(), nperm = 1000, mrank = FALSE)
{
# MRM: Multiple Regression on distance Matrices
# Sarah Goslee 2008-07-18
# tests R^2 and regression coefficients using a
# permutation test

# Stuff R needs to be able to use a formula
     m <- match.call(expand.dots = FALSE)
     m2 <- match(c("formula", "data"), names(m), nomatch=0)
     m <- m[c(1, m2)]
     m[[1]] <- as.name("model.frame")
     m <- eval(m, parent.frame())
     m <- as.matrix(m)


# End of R stuff. m is now the data for the MRM test as
# columns y, x, n1, n2, n3, ...
# Determine the size of the matrices & do some error checking.
     n <- (1 + sqrt(1 + 8 * nrow(m)))/2
     if(abs(n - round(n)) > 0.0000001)
stop("Matrix not square.\n")
     n <- round(n)
     if(ncol(m) < 2) stop("Not enough data. \n")

    if(mrank) {
        m <- apply(m, 2, rank)
    }

    # convert matrices to column order to ensure compatibility with C
    for(thiscol in 1:ncol(m)) {
        tempmat <- full(m[,thiscol])
        m[,thiscol] <- tempmat[col(tempmat) > row(tempmat)]
    }


    # use matrix form to speed up calculations
    X <- m[ ,2:ncol(m), drop=FALSE]
    X <- cbind(rep(1, nrow(X)), X)
    Y <- m[ ,1, drop=FALSE]

    nd <- nrow(X)

    # only need to calculate (X'X)^-1 once
    XX <- crossprod(X)
    XX <- solve(XX)

    # will need to calculate Xy for each permutation
    XY <- crossprod(X, Y)
    YY <- crossprod(Y)

    # regression coefficients
    b <- XX %*% XY
    rownames(b) <- c("Int", colnames(X)[2:ncol(X)])

    bXY <- crossprod(b, XY)
    SSE <- YY - bXY

    SSTO <- YY - sum(Y)^2/nd
    SSR = SSTO - SSE

    # R2 = 1 - SSE/SSTO
    R2 <- 1 - SSE/SSTO
    R2 <- as.vector(R2)

    # F* = MSR / MSE
    # MSR = SSR / (p - 1) 
    # MSE = SSE / (n - p)
    p <- ncol(X) # number of parameters estimated
    F <- (SSR / (p - 1)) / (SSE / (nd - p))

    R2.pval <- NA
    b.pval <- rep(NA, ncol(X))
    F.pval <- NA

    if(nperm > 0) {

        R2.all <- numeric(nperm)

        # for regression coefficients, use pseudo-t of Legendre et al. 1994
        b.all <- numeric(nperm*p)

        # try out an overall F-test for lack of fit
        F.all <- numeric(nperm)

        cresults <- .C("mrmperm", as.double(as.vector(X)), as.double(as.vector(Y)), 
             as.integer(p), as.integer(nd), as.integer(n), as.integer(nperm), 
             R2.all = as.double(R2.all), b.all = as.double(b.all), 
             F.all = as.double(F.all), as.double(numeric(n*n)), 
             as.integer(numeric(n)), as.double(as.vector(XX)), as.double(numeric(p)),
             as.double(0), as.double(numeric(p)), PACKAGE = "ecodist")

        R2.all <- cresults$R2.all
        R2.pval <- length(R2.all[R2.all >= R2.all[1]])/nperm

        F.all <- cresults$F.all
        F.pval <- length(F.all[F.all >= F.all[1]])/nperm

        # b.all contains pseudo-t of Legendre et al. 1994
        b.all <- matrix(cresults$b.all, nrow=nperm, ncol=p, byrow=TRUE)
        b.pval <- apply(b.all, 2, function(x)length(x[abs(x) >= abs(x[1])])/nperm)
}

    list(coef=cbind(b, pval=b.pval), r.squared=c(R2=R2, pval = R2.pval),F.test=c(F=F, F.pval = F.pval))
}

