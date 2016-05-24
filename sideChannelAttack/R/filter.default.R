cor2I2 <-
function (rho) 
{
    rho <- pmin(rho, 1 - 1e-05)
    -1/2 * log(1 - rho^2)
}
filter.NULL.default <-
function (...) 
{
    res = list()
    class(res) <- "filter.NULL"
    return(res)
}
filter.PCA.default <-
function (X, nbreVarX_, ...) 
{
    if (!is.matrix(X)) {
        stop("'X' has to be a matrix")
        return(-1)
    }
    if (nbreVarX_ > dim(X)[2]) {
        stop("the number of variable to take from 'X' has to be less than the number of variable of 'X'")
        return(-1)
    }
    if (nbreVarX_ <= 0) {
        stop("the number of variable to take from 'X' has to be positive")
        return(-1)
    }
    res = princomp(x=X)
    res2 = list(mod = res, nbreVarX = nbreVarX_)
    class(res2) <- "filter.PCA"
    return(res2)
}
filter.RegressionTreeFilter.default <-
function (X, nbreVarX_, ...) 
{
    if (nbreVarX_ <= 0) {
        stop("the number of variable to take from 'X' has to be positive")
        return(-1)
    }
    if (nbreVarX_ > dim(X)[2]) {
        stop("the number of variable to take from 'X' has to be less than the number of variable of 'X'")
        return(-1)
    }
    res = list(nbreVarX = nbreVarX_)
    class(res) <- "filter.RegressionTreeFilter"
    return(res)
}
filter.mRMR.default <-
function (X, Y, nbreVarX_, ...) 
{
    if (!is.matrix(X)) {
        stop("'X' has to be a matrix")
        return(-1)
    }
    if (!is.vector(Y)) {
        stop("'Y' has to be a vector")
        return(-1)
    }
    if (nbreVarX_ > dim(X)[2]) {
        stop("the number of variable to take from 'X' has to be less than the number of variable of 'X'")
        return(-1)
    }
    if (nbreVarX_ <= 0) {
        stop("the number of variable to take from 'X' has to be positive")
        return(-1)
    }
    if (dim(X)[1] != length(Y)) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
	mim <- mutinformation(cbind(discretize(X),Y))
	n <- dim(X)[2]
	S <- rep(0,n)
	Best <- c()
	TailleS <- 0
	for(i in 1:n) {
		iIndiceBest <- -1
		iCoutBest <- -1
		for(j in which(S == 0)) {
			jCout <- mim[n+1,j] - ifelse(TailleS == 0,0,mean(mim[j,which(S==1)]))
			if(jCout > iCoutBest){
				iCoutBest <- jCout
				iIndiceBest <- j
			}
		}
		Best <- c(Best,iIndiceBest)
		S[iIndiceBest] <- 1
		TailleS <- TailleS + 1
		if(TailleS >= nbreVarX_)
			break;
	}
	res <- list(filter = Best)
    class(res) <- "filter.mRMR"
    return(res)
}


filter.MAX.default <-
function (nbreVarX_,...) 
{
    res = list(nbreVarX = nbreVarX_)
    class(res) <- "filter.MAX"
    return(res)
}
