# genData.R -- version 2010-12-29
genData <- function(nP, nO, ol, dy) {
    # create data file [Salibian-Barrera & Yohai 2006]
    # nP = 36  % regressors
    # nO = 400 % number of obs
    # ol = 40  % number of outliers 
    # dy = 110 % outlier size ('M' in S-B&Y 2006): 90 to 200
    mRN <- function(m, n) array(rnorm(m * n), dim = c(m, n))
    y <- mRN(nO,1L)
    X <- cbind(as.matrix(numeric(nO) + 1), mRN(nO, nP - 1L))
    zz <- sample(nO)
    z <- cbind(1,100, array(0, dim = c(1L, nP - 2L)))
    for (i in 1L:ol){
        X[zz[i], ] <- z
        y[zz[i]] <- dy
    }
    list(X = X, y = y)
}