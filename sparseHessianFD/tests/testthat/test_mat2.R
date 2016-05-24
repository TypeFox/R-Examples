## Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun



context("Matrix.to.Pointers")
test_that("Matrix.to.Pointers", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,1)
    cols <- c(1,2,2,3,4,4,5)
    vals <- seq(10, 10+nnz-1)
    A <- sparseMatrix(i=rows, j=cols, x=vals, dims=c(k,k), index1=TRUE)
    AR <- as(A, "RsparseMatrix")
    AC <- as(A, "CsparseMatrix")
    AT <- as(A, "TsparseMatrix")

    P1r <- Matrix.to.Pointers(A, FALSE, TRUE, "row", TRUE)
    P1c <- Matrix.to.Pointers(A, FALSE, TRUE, "column", TRUE)
    P1t <- Matrix.to.Pointers(A, FALSE, TRUE, "triplet", TRUE)

    expect_equal(names(P1r), c("jCol","ipntr","x","class"))
    expect_equal(names(P1c), c("iRow","jpntr","x", "class"))
    expect_equal(names(P1t), c("rows","cols","x","class"))

    A1r <- as(sparseMatrix(j=P1r$jCol, p=P1r$ipntr-1, x=P1r$x),"RsparseMatrix")
    A1c <- sparseMatrix(i=P1c$iRow, p=P1c$jpntr-1, x=P1c$x)
    A1t <- sparseMatrix(i=P1t$rows, j=P1t$cols, x=P1t$x, giveCsparse=FALSE)

    expect_equal(AR, A1r)
    expect_equal(AC, A1c)
    expect_equal(uniqTsparse(AT), uniqTsparse(A1t))

    P2r <- Matrix.to.Pointers(AR, FALSE, TRUE, "row", TRUE)
    P2c <- Matrix.to.Pointers(AR, FALSE, TRUE, "column", TRUE)
    P2t <- Matrix.to.Pointers(AR, FALSE, TRUE, "triplet", TRUE)

    expect_equal(names(P2r), c("jCol","ipntr","x","class"))
    expect_equal(names(P2c), c("iRow","jpntr","x", "class"))
    expect_equal(names(P2t), c("rows","cols","x","class"))

    A2r <- as(sparseMatrix(j=P2r$jCol, p=P2r$ipntr-1, x=P2r$x),"RsparseMatrix")
    A2c <- sparseMatrix(i=P2c$iRow, p=P2c$jpntr-1, x=P2c$x)
    A2t <- sparseMatrix(i=P2t$rows, j=P2t$cols, x=P2t$x, giveCsparse=FALSE)

    expect_equal(AR, A2r)
    expect_equal(AC, A2c)
    expect_equal(uniqTsparse(AT), uniqTsparse(A2t))

    P3r <- Matrix.to.Pointers(AC, FALSE, FALSE, "row", TRUE)
    P3c <- Matrix.to.Pointers(AC, FALSE, FALSE, "column", TRUE)
    P3t <- Matrix.to.Pointers(AC, FALSE, FALSE, "triplet", TRUE)

    expect_equal(names(P3r), c("jCol","ipntr","class"))
    expect_equal(names(P3c), c("iRow","jpntr", "class"))
    expect_equal(names(P3t), c("rows","cols","class"))

    A3r <- as(sparseMatrix(j=P3r$jCol, p=P3r$ipntr-1),"RsparseMatrix")
    A3c <- sparseMatrix(i=P3c$iRow, p=P3c$jpntr-1)
    A3t <- sparseMatrix(i=P3t$rows, j=P3t$cols, giveCsparse=FALSE)

    nAR <- as(as(AR,"nMatrix"),"RsparseMatrix")
    nAC <- as(as(AC,"nMatrix"),"CsparseMatrix")
    nAT <- as(as(AT,"nMatrix"),"TsparseMatrix")

    expect_equal(nAR, A3r)
    expect_equal(nAC, A3c)
    expect_equal(uniqTsparse(nAT), uniqTsparse(A3t))

})
