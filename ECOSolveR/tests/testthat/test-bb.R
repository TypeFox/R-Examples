context("ECOS BB Tests")
## test 1a_bool
test_that("test1a returns the right answer", {
    G <- local({
        Gpr <- c(2.0, 3.0, 1.0, 4.0)
        Gjc <- as.integer(c(0, 2, 4))
        Gir <- as.integer(c(0, 1, 0, 1))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })
    c <- c(-1.1, -1.)
    h <- c(4., 12.)
    dims <- list(l=2L, q=NULL, e=0L)
    bool_idx <- 1L

    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)

    ## Answer:
    ## pfloat x[2] = {1.0, 2.0};
    expect_equal(c(1.0, 2.0), retval$x)
})


## test 1a_int
test_that("test1a_int returns the right answer", {
    G <- local({
        Gpr <- c(2.0, 3.0, -1, 1.0, 4.0, -1)
        Gjc <- as.integer(c(0, 3, 6))
        Gir <- as.integer(c(0, 1, 2, 0, 1, 3))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })

    c <- c(-1., -1.1)
    h <- c(4., 12., 0. , 0.)
    dims <- list(l=4L, q=NULL, e=0L)
    int_idx <- as.integer(c(1,2))

    retval <- ECOS_csolve(c = c, G = G, h = h,
                          dims = dims,
                           A = NULL, b = numeric(0),
                           int_vars = int_idx)
    ## Answer:
    ## pfloat x[2] = {0.0, 3.0};
    expect_equal(c(0.0, 3.0), retval$x)
})


## test 1b
test_that("test1b returns the right answer", {
    G <- local({
        Gpr <- c(2.0, 3.0, 1.0, 4.0)
        Gjc <- as.integer(c(0, 2, 4))
        Gir <- as.integer(c(0, 1, 0, 1))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })

    c <- c(-1.0, -1.0)
    h <- c(4.0, 12.0)
    dims <- list(l=2L, q=NULL, e=0L)
    bool_idx <- 2L
    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)
    ## Answer:
    ## pfloat x[2] = {1.5, 1.0};
    expect_equal(c(1.5, 1.0), retval$x)
})

## test 2
test_that("test2 returns the right answer", {
    G <- local({
        Gpr <- c(2.0, 3.0, 1.0, 4.0)
        Gjc <- as.integer(c(0, 2, 4))
        Gir <- as.integer(c(0, 1, 0, 1))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })
    c <- c(-1.0, -1.0)
    h <- c(4.0, 12.0)
    dims <- list(l=2L, q=NULL, e=0L)
    bool_idx <- as.integer(c(1,2))
    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)
    ## Answer:
    ## pfloat x[2] = {1.0, 1.0};
    expect_equal(c(1.0, 1.0), retval$x)
})

## test 3
test_that("test3 returns the right answer", {
    G <- local({
        Gpr <- c(2.0, 3.0, 1.0, 4.0)
        Gjc <- as.integer(c(0, 2, 4))
        Gir <- as.integer(c(0, 1, 0, 1))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })
    c <- c(1.0, -1.0)
    h <- c(4.0, 12.0)
    dims <- list(l=2L, q=NULL, e=0L)
    bool_idx <- 1L
    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)
    ## Answer:
    ## pfloat x[2] = {0.0, 3.0};
    expect_equal(c(0.0, 3.0), retval$x)
})


## test 4
test_that("test4 returns the right answer", {
    G <- local({
        Gpr <- c(2.0,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1)
        Gjc <- as.integer(c(0,3,6,9,12,15,18))
        Gir <- as.integer(c(0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })
    c <- c(3.0, 5, 6, 9, 10, 10)
    h <- c(-2.0, 2, -3)
    dims <- list(l=3L, q=NULL, e=0L)
    bool_idx <- as.integer(c(0,1,2,3,4,5) + 1)
    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)
    ## Answer:
    ## pfloat x[6] = {0,1,1,0,0,0};
    expect_equal(c(0.0, 1.0, 1.0, 0.0, 0.0, 0.0), retval$x)
})

## test 5
test_that("test4 returns the right answer", {
    G <- local({
        Gpr <- c(2.0,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1)
        Gjc <- as.integer(c(0,3,6,9,12,15,18))
        Gir <- as.integer(c(0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2))
        sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
    })
    c <- c(3.0, 5, 6, 9, -10, 10)
    h <- c(-2.0, 2, -3)
    dims <- list(l=3L, q=NULL, e=0L)
    bool_idx <- as.integer(c(0,1,2,3,4,5) + 1)
    retval <- ECOS_csolve(c = c, G = G, h = h,
                           dims = dims,
                           A = NULL, b = numeric(0),
                           bool_vars = bool_idx)
    ## Answer:
    ## pfloat x[6] = {0,0,1,1,1,0};
    expect_equal(c(0.0, 0.0, 1.0, 1.0, 1.0, 0.0), retval$x)
})


##
## Ignore tests 6 and 7 since 6 is supposed to fail and it is not
## clear test 7 has been vetted carefully...
##
## test 6
## G <- local({
##     Gpr <- c(2.0,5,-5,-6,3,1,3,-1,-4,-4,-3,2,-1,2,-2,2,-1,1)
##     Gjc <- as.integer(c(0,3,6,9,12,15,18))
##     Gir <- as.integer(c(0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2))
##     sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
## })
## c <- c(-3, -5, 6, 9, 10, -10.0)
## h <- c(-2.0, 1, -3)
## dims <- list(l=3L, q=NULL, e=0L)
## bool_idx <- as.integer(c(0,1,2,3,4,5) + 1)
## retval <- ECOS_csolve(c = c, G = G, h = h,
##                        dims = dims,
##                        A = NULL, b = numeric(0),
##                        bool_vars = bool_idx)

## sprintf("%f", retval$x)

## ## test 7

## G <- local({
##     Gpr <- c(9999,-9999,9999,-9999,9999,-9999,9999,-9999,9999,-9999,-3.5008,3.5008,-0.4504,0.4504,-0.8764999999999999,0.8764999999999999,-0.1088,0.1088,1,1,-1,-8.4095,8.4095,-1.0107,1.0107,-1.686,1.686,-0.3525,0.3525,1,1,-1,-15.1987,15.1987,-2.0203,2.0203,-2.3932,2.3932,-0.6233,0.6233,1,1,-1,-22.5405,22.5405,-3.1862,3.1862,-2.8749,2.8749,-0.7923,0.7923,1,1,-1,-29.2639,29.2639,-4.3096,4.3096,-3.0189,3.0189,-0.8116,0.8116,1,1,-1,3.5008,-3.5008,0.4504,-0.4504,0.8764999999999999,-0.8764999999999999,0.1088,-0.1088,1,1,-1,8.4095,-8.4095,1.0107,-1.0107,1.686,-1.686,0.3525,-0.3525,1,1,-1,15.1987,-15.1987,2.0203,-2.0203,2.3932,-2.3932,0.6233,-0.6233,1,1,-1,22.5405,-22.5405,3.1862,-3.1862,2.8749,-2.8749,0.7923,-0.7923,1,1,-1,29.2639,-29.2639,4.3096,-4.3096,3.0189,-3.0189,0.8116,-0.8116,1,1,-1)
##     Gjc <- as.integer(c(0,2,4,6,8,10,21,32,43,54,65,76,87,98,109,120))
##     Gir <- as.integer(c(8,9,10,11,12,13,14,15,16,17,0,1,2,3,4,5,6,7,8,18,19,0,1,2,3,4,5,6,7,10,18,20,0,1,2,3,4,5,6,7,12,18,21,0,1,2,3,4,5,6,7,14,18,22,0,1,2,3,4,5,6,7,16,18,23,0,1,2,3,4,5,6,7,9,18,24,0,1,2,3,4,5,6,7,11,18,25,0,1,2,3,4,5,6,7,13,18,26,0,1,2,3,4,5,6,7,15,18,27,0,1,2,3,4,5,6,7,17,18,28))
##     sparseMatrix(x=Gpr, i=Gir, p=Gjc, index1=FALSE)
## })
## c <- c(0,0,0,0,0,0.127,0.9134,0.6324,0.0975,0.2785,0.873,0.0866,0.3676,0.9025,0.7215)
## h <- c(-729.9349999999999,789.9349999999999,-71.015,131.015,-89.66,149.66,-1.165,61.165,9999,0,9999,0,9999,0,9999,0,9999,0,150,0,0,0,0,0,0,0,0,0,0)
## dims <- list(l=29L, q=NULL, e=0L)
## bool_idx <- as.integer(c(0, 1, 2, 3, 4) + 1)
## retval <- ECOS_csolve(c = c, G = G, h = h,
##                        dims = dims,
##                        A = NULL, b = numeric(0),
##                        bool_vars = bool_idx)

## sprintf("%f", retval$x)
