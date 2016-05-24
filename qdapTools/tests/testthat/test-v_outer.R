context("Checking v_outer")

test_that("v_outer works on data.frames",{
    #|------------------------------------------------------|
    #|    SETTING UP VARIOUS FUNCTIONS THAT WILL BE USED    |
    #|------------------------------------------------------|
    pooled.sd <- function(x, y) {
        n1 <- length(x)
        n2 <- length(y)
        s1 <- sd(x)
        s2 <- sd(y)
        sqrt(((n1-1)*s1 + (n2-1)*s2)/((n1-1) + (n2-1)))
    }
    
    euc.dist <- function(x,y) sqrt(sum((x - y) ^ 2))
    sum2 <- function(x, y) sum(x, y)
    arbitrary <- function(x, y) round(sqrt(sum(x)) - sum(y), digits=1)
	
    ## Cosine similarity
    cos_sim <- function(x, y) x %*% y / sqrt(x%*%x * y%*%y)
    #--------------------------------------------------------#

    v <- v_outer(mtcars[, 1:4], cor)
    w <- v_outer(mtcars[, 1:4], pooled.sd)
    x <- v_outer(mtcars[, 1:4], euc.dist)
    y <- v_outer(mtcars[, 1:4], sum2)
    z <- v_outer(mtcars[, 1:4], arbitrary)
    
    expect_true(isSymmetric(w))
    expect_false(isSymmetric(z))
    L1 <- list(v, w, x, y, z)
    expect_false(any(sapply(L1, is.data.frame)))
    expect_true(all(sapply(L1, is.matrix)))
    expect_true(sum(v - cor(mtcars[, 1:4])) == 0)
	

})


test_that("v_outer works on lists",{
    #|------------------------------------------------------|
    #|    SETTING UP VARIOUS FUNCTIONS THAT WILL BE USED    |
    #|------------------------------------------------------|
    pooled.sd <- function(x, y) {
        n1 <- length(x)
        n2 <- length(y)
        s1 <- sd(x)
        s2 <- sd(y)
        sqrt(((n1-1)*s1 + (n2-1)*s2)/((n1-1) + (n2-1)))
    }
    
    euc.dist <- function(x,y) sqrt(sum((x - y) ^ 2))
    sum2 <- function(x, y) sum(x, y)
    arbitrary <- function(x, y) round(sqrt(sum(x)) - sum(y), digits=1)
	    
    ## Cosine similarity
    cos_sim <- function(x, y) x %*% y / sqrt(x%*%x * y%*%y)
    #--------------------------------------------------------#

    mtcars2 <- lapply(mtcars[, 1:4], "[")
    v <- v_outer(mtcars2, cor)
    w <- v_outer(mtcars2, pooled.sd)
    x <- v_outer(mtcars2, euc.dist)
    y <- v_outer(mtcars2, sum2)
    z <- v_outer(mtcars2, arbitrary)
    
    expect_true(isSymmetric(w))
    expect_false(isSymmetric(z))
    L1 <- list(v, w, x, y, z)
    expect_false(any(sapply(L1, is.data.frame)))
    expect_true(all(sapply(L1, is.matrix)))
    expect_true(sum(v - cor(mtcars[, 1:4])) == 0)

})


test_that("v_outer works on matrices",{
    #|------------------------------------------------------|
    #|    SETTING UP VARIOUS FUNCTIONS THAT WILL BE USED    |
    #|------------------------------------------------------|
    pooled.sd <- function(x, y) {
        n1 <- length(x)
        n2 <- length(y)
        s1 <- sd(x)
        s2 <- sd(y)
        sqrt(((n1-1)*s1 + (n2-1)*s2)/((n1-1) + (n2-1)))
    }
    
    euc.dist <- function(x,y) sqrt(sum((x - y) ^ 2))
    sum2 <- function(x, y) sum(x, y)
    arbitrary <- function(x, y) round(sqrt(sum(x)) - sum(y), digits=1)
	    
    ## Cosine similarity
    cos_sim <- function(x, y) x %*% y / sqrt(x%*%x * y%*%y)
    #--------------------------------------------------------#

    set.seed(10)
    mat <- matrix(rbinom(50, 0:1, .45), ncol=5)
    v <- v_outer(mat, cor)
    w <- v_outer(mat, pooled.sd)
    x <- v_outer(mat, euc.dist)
    y <- v_outer(mat, sum2)
    z <- v_outer(mat, arbitrary)
    
    expect_true(isSymmetric(w))
    expect_false(isSymmetric(z))
    L1 <- list(v, w, x, y, z)
    expect_false(any(sapply(L1, is.data.frame)))
    expect_true(all(sapply(L1, is.matrix)))
    expect_true(sum(v - cor(mat)) == 0)
	

})