context("corbetw2mat")

test_that('corbetw2mat works with what="paired"', {

    n_ind <- 20
    n_col <- 3

    x <- matrix(rnorm(n_ind*n_col), ncol=n_col)
    y <- x + matrix(rnorm(n_ind*n_col, 0, 0.5), ncol=n_col)
    colnames(x) <- colnames(y) <- 1:n_col
    rownames(x) <- rownames(y) <- 1:n_ind

    result <- corbetw2mat(x, y, what="paired")

    expected <- NULL
    for(i in 1:ncol(x))
        expected[i] <- cor(x[,i], y[,i])
    names(expected) <- 1:n_col

    expect_equal(result, expected)

})

test_that('corbetw2mat works in the other cases', {

    n_ind <- 20
    n_col <- 3
    n_col_extra <- 5

    set.seed(7490138)

    x <- matrix(rnorm(n_ind*n_col), ncol=n_col)
    y <- cbind(x + rnorm(n_ind*n_col, 0, 0.5),
               matrix(rnorm(n_ind*n_col_extra, 0, 0.5), ncol=n_col_extra))
    y <- y[,sample(ncol(y))]
    colnames(x) <- 1:ncol(x)
    colnames(y) <- 1:ncol(y)
    rownames(x) <- rownames(y) <- 1:n_ind

    result1 <- corbetw2mat(x, y, what="bestright")
    expected1 <- data.frame(cor=c(0.9360757383422467,0.9513932869294095,0.7754314786368088),
                            yindex=c(7,2,6),
                            ycol=c("7","2","6"))
    rownames(expected1) <- colnames(x)
    expect_equal(result1, expected1)


    result2 <- corbetw2mat(x, y, what="bestpairs", corthresh = 0.3)
    expected2 <- data.frame(cor=c(0.3779868370283727,0.9360757383422467,0.9513932869294095,0.7754314786368088),
                            xindex=c(1,1,2,3),
                            yindex=c(4,7,2,6),
                            xcol=c("1","1","2","3"),
                            ycol=c("4","7","2","6"))
    rownames(expected2) <- 1:nrow(expected2)
    expect_equal(result2, expected2)

    result3 <- corbetw2mat(x, y, what="all")
    expected3 <- cor(cbind(x,y))[1:ncol(x), ncol(x) + 1:ncol(y)]
    expect_equal(result3, expected3)

})
