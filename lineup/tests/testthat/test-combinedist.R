context("combinedist")

test_that('combinedist works in a simple case', {

    # simulate MVN, 100 individuals, 40 measurements (of which 20 are just noise)
    V <- matrix(0.3, ncol=20, nrow=20) + diag(rep(0.5, 20))
    D <- chol(V)
    z <- matrix(rnorm(20*100), ncol=20) %*% D

    # create three data matrices as z + noise
    x <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
    y <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))
    w <- cbind(z + rnorm(20*100, 0, 0.2), matrix(rnorm(20*100), ncol=20))

    # permute some rows of x
    x[51:53,] <- x[c(52,53,51),]

    # add column and row names
    dimnames(x) <- dimnames(y) <- dimnames(w) <-
        list(paste("ind", 1:100, sep=""), paste("gene", 1:40, sep=""))

    # calculate correlations between cols of x and of the other two matrices
    corxy <- corbetw2mat(x, y)
    corxw <- corbetw2mat(x, w)

    # using columns with corr > 0.75,
    # calculate distance (using "correlation" as a measure...really similarity)
    dxy <- distee(x[,corxy>0.75], y[,corxy>0.75], d.method="cor", labels=c("x", "y"))
    dxw <- distee(x[,corxw>0.75], w[,corxw>0.75], d.method="cor", labels=c("x", "w"))

    expect_equivalent(combinedist(dxy, dxw), (dxy+dxw)/2)

})
