## Part of the sparseHessianFD package
## Copyright (C) 2013-2015 Michael Braun

context("Matrix.to.Coord")
test_that("Matrix.to.Coord", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,1)
    cols <- c(1,2,2,3,4,4,5)

    M1.true <- sparseMatrix(i=rows, j=cols, dims=c(k,k))
    C1 <- Matrix.to.Coord(M1.true)
    expect_equal(names(C1), c("rows","cols"))
    M1 <- sparseMatrix(i=C1$rows, j=C1$cols)
    expect_equal(M1.true, M1)
})


context("Matrix.to.Pointers")
test_that("Matrix.to.Pointers", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,1)
    cols <- c(1,2,2,3,4,4,5)

    M1.true <- sparseMatrix(i=rows, j=cols, dims=c(k,k) )
    C1r <- Matrix.to.Pointers(M1.true, order="row", index1=TRUE)
    C1c <- Matrix.to.Pointers(M1.true, order="column", index1=TRUE)

    expect_equal(names(C1r), c("jCol","ipntr","class"))
    expect_equal(names(C1c), c("iRow","jpntr","class"))

    M1r <- sparseMatrix(j=C1r$jCol, p=C1r$ipntr-1, dims=c(k,k))
    M1c <- sparseMatrix(i=C1c$iRow, p=C1c$jpntr-1, dims=c(k,k))
    expect_equal(M1.true, M1r)
    expect_equal(M1.true, M1c)

    M2.true <- sparseMatrix(i=rows, j=cols, dims=c(k,k))
    C2r <- Matrix.to.Pointers(M2.true, order="row", index1=FALSE)
    C2c <- Matrix.to.Pointers(M2.true, order="column", index1=FALSE)
    M2r <- sparseMatrix(j=C2r$jCol, p=C2r$ipntr, dims=c(k,k), index1=FALSE)
    M2c <- sparseMatrix(i=C2c$iRow, p=C2c$jpntr, dims=c(k,k), index1=FALSE)

    expect_equal(M2.true, M2r)
    expect_equal(M2.true, M2c)
})

context("Coord.to.Pointers")
test_that("Coord.to.Pointers", {

    ## LT matrix data
    nnz <- 7
    k <- 5
    rows <- c(1,2,5,3,4,5,2)
    cols <- c(1,2,2,3,4,5,5)

    M1.true <- sparseMatrix(i=rows, j=cols, dims=c(k,k) )
    C1s <- Coord.to.Pointers(rows, cols, c(k,k), triangle=FALSE,
                             symmetric=TRUE,
                             order="column")
    C1r <- Coord.to.Pointers(rows, cols, c(k,k), triangle=FALSE,
                             order="row")
    C1c <- Coord.to.Pointers(rows, cols, c(k,k), triangle=FALSE,
                             order="column")



    M1s <- sparseMatrix(i=C1s$idx, p=C1s$pntr-1)
    M1r <- sparseMatrix(j=C1r$jCol, p=C1r$ipntr-1)
    M1c <- sparseMatrix(i=C1c$iRow, p=C1c$jpntr-1)

    expect_equal(M1.true, M1s)
    expect_equal(M1.true, M1r)
    expect_equal(M1.true, M1c)

    rowsLT <- c(1,2,5,3,4,4,5)
    colsLT <- c(1,2,2,3,3,4,5)

    M2.true <- as(tril(sparseMatrix(i=rowsLT, j=colsLT, dims=c(k,k))),"ngCMatrix")
    M2.trueS <- as(sparseMatrix(i=rowsLT, j=colsLT, dims=c(k,k), symmetric=TRUE),"ngCMatrix")
    expect_true(Matrix::isTriangular(M2.true))
    C2s <- Coord.to.Pointers(rowsLT, colsLT, c(k,k), triangle=TRUE,
                             symmetric=TRUE, order="column")
    C2r <- Coord.to.Pointers(rowsLT, colsLT, c(k,k), triangle=TRUE,
                             symmetric=FALSE, order="row")
    C2c <- Coord.to.Pointers(rowsLT, colsLT, c(k,k), triangle=TRUE,
                             symmetric=FALSE, order="column")

    M2s <- sparseMatrix(i=C2s$idx, p=C2s$pntr-1)
    M2r <- sparseMatrix(j=C2r$jCol, p=C2r$ipntr-1)
    M2c <- sparseMatrix(i=C2c$iRow, p=C2c$jpntr-1)

    expect_true(Matrix::isSymmetric(M2s))
    expect_equal(M2.trueS, M2s)
    expect_equal(M2.true, M2r)
    expect_equal(M2.true, M2c)

    P1 <- Coord.to.Pattern.Matrix(rows, cols, dims=c(k,k),
                                  symmetric=FALSE)

    P2 <- Coord.to.Pattern.Matrix(rowsLT, colsLT, dims=c(k,k),
                                  symmetric=FALSE)

    P3 <- as(Coord.to.Pattern.Matrix(rowsLT, colsLT, dims=c(k,k),
                                  symmetric=TRUE), "ngCMatrix")

    expect_equal(M1.true, P1)
    expect_equal(M2.true, P2)
    expect_equal(M2.trueS, P3)
})
