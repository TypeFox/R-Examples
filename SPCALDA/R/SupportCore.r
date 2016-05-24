


mypreparedata = function (X, Y, rescale = FALSE)
{
#basic check
    if (is.null(dim(X)))
        stop("\n'variables' must be a matrix or data.frame")
    if (nrow(X) != length(Y))
        stop("\n'variables' and 'labels' must have same lengths")
# label as factors
    if (!is.vector(Y) && !is.factor(Y))
        stop("\n'group' must be a factor")
    if (!is.factor(Y))
        Y = as.factor(Y)
#column and row names
    if (is.null(colnames(X)))
        colnames(X) = paste(rep("X", ncol(X)), seq_len(ncol(X)),
            sep = "")
    if (is.null(rownames(X)))
        rownames(X) = 1L:nrow(X)
#scale data when necessary (temply, let us not rescale the data)
    if (rescale) X=scale(X)

    return(list(X = X, Y = Y))
}



# caluclate NewX which is a (n+g)*p matrix; ##it can replace AllCov in the future
NewX = function(X,Y){
    n = nrow(X); p = ncol(X)
    if (!is.factor(Y)) Y = as.factor(Y)
    Glevels = levels(Y)
    G = nlevels(Y)

    AllMean = colMeans(X)
    NewX2 = NULL
    for (g in 1:G) {
        index = which(Y == Glevels[g])
        n_g = length(index)
        mu_g = colMeans(X[index, ])
        X[index,] = X[index,] - matrix(mu_g,n_g,p,byrow=TRUE)
        NewX2= rbind(NewX2,sqrt(n_g)*(mu_g-AllMean))
    }
    return(list(NewX1=X,NewX2=NewX2))
}

RRRotation = function(X,Y,rho,K=min(dim(X))){   # input X=n\times p data matrix, Y label, rho = values of tuning parameters, K = number of output PC
    ob = mypreparedata(X, Y, rescale = FALSE)
    X = ob$X; Y = ob$Y
    n = nrow(X); p = ncol(X); M = length(rho)
    RotationMatrix = list()
    if (p<n) { # compute W + rho B
       object = NewX(X,Y)
       X1 = object$NewX1; X2 = object$NewX2
       W = crossprod(X1)/n; B = crossprod(X2)/n
       for (m in 1:M) {
          Sigma_rho = W + rho[m]*B
          # SVD or eigen value decomposition
          currentRotationMatrix = (eigen(Sigma_rho, symmetric=TRUE)$vectors)[,1:K]
          RotationMatrix = c(RotationMatrix,list(currentRotationMatrix))
       }
       return(RotationMatrix)
    }
    else {
       object = NewX(X,Y)
       X1 = object$NewX1; X2 = object$NewX2
       for (m in 1:M) {
          NewXm = rbind(X1,sqrt(rho[m])*X2)
          small_Sigma_rho = tcrossprod(NewXm)/n
          # SVD or eigen value decomposition
          temp = (eigen(small_Sigma_rho, symmetric=TRUE)$vectors)[,1:K]
          currentRotationMatrix = scale(t(NewXm)%*%temp,center=FALSE)/sqrt(p-1)
          RotationMatrix = c(RotationMatrix,list(currentRotationMatrix))
       }
       return(RotationMatrix)
    }
}



 

nestedLDA = function(x,y,xt,yt) {
     Errors = NULL
     for (i in 1:dim(x)[2]) {
        if (i==1) {
           x1=as.matrix(x[,1]); x2 = as.matrix(xt[,1])}
        else
           {   x1=x[,1:i]; x2 = xt[,1:i]}
        ob= lda(x1,y)
        yhat = predict(ob, x2)$class
        er = sum(yhat != yt)
        Errors = c(Errors, er)
     }
     return(Errors)
}

# return the position of one of minimal element in a matrix
tiebreaker = function(A, smallestK=TRUE){
	if (length(which(A == min(A)))==1) return(arrayInd(which.min(A), dim(A)))
	ob = which(A == min(A), arr.ind = TRUE) # these are positions of all minimal elements
	# find the ROW (upper/smaller one if tie) that includes most minimal elements
	colnum = ob[,1]
	Ai = which.max(tabulate(colnum))
	# in that ROW, find the position of minimal element (left/smaller one if tie)
    if (smallestK) {Aj = which.min(A[Ai,])}
    # alternative way to determine Aj (median position)
    else {temp = min(A[Ai,]); Aj = floor(median(which(A[Ai,]==temp)))}
    return(c(Ai,Aj))
}
