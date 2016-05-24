predict.sPLS <- predict.gPLS <- predict.sgPLS <- function(object, newdata,  ...)
{
    
    #-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
    
    X = object$X
    Y = object$Y
    q = ncol(Y)
    p = ncol(X)
    
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
        " or a vector of length = ", p, ".")
    }
    
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
        " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p)
    }
    
    #-- initialisation des matrices --#
    ncomp = object$ncomp
    a = object$loadings$X
    b = object$loadings$Y
    c = object$mat.c
    
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
    
    newdata = as.matrix(newdata)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    Y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
    
    variates.X = object$variates$X
    betay = list()
    
    #-- prediction --#
    for(h in 1:ncomp){
        
        dd= coefficients(lm(Y~variates.X[,1:h,drop=FALSE])) #regression of Y on variates.global.X => =loadings.global.Y at a scale factor
        if(q==1){betay[[h]]=(dd[-1])}
        if(q>=2){betay[[h]]=(dd[-1,])}
        
        W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
        B = W %*% drop(betay[[h]])
        
        Y.temp=scale(newdata,center=means.X,scale=sigma.X) %*% as.matrix(B) #so far: gives a prediction of Y centered and scaled
        Y.temp2=scale(Y.temp,center=FALSE,scale=1/sigma.Y) #so far: gives a prediction of Y centered, with the right scaling
        Y.temp3=scale(Y.temp2,center=-means.Y,scale=FALSE) #so far: gives a prediction of Y with the right centering and scaling
        
        Y.hat[, , h] = Y.temp3 # we add the variance and the mean of Y used in object to predict
        t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
        B.hat[, , h] = B
    }  #end h
    
    #-- valeurs sortantes --#
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
    
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat,betay=betay)))
}


