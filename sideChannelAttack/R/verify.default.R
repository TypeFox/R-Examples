verify.loo.default <-
function (model, filter, X, Y, nbreVarX, param.model=list(), param.fs=list(), ...) 
{
    if (!is.vector(Y)) {
        stop("'Y' has to be a vector")
        return(-1)
    }
    if (!is.matrix(X)) {
        stop("'X' has to be a matrix")
        return(-1)
    }
    if (length(Y) != dim(X)[1]) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
    TP <- rep(0, length(nbreVarX))
    TN <- rep(0, length(nbreVarX))
    FN <- rep(0, length(nbreVarX))
    FP <- rep(0, length(nbreVarX))
    nbreVarX_dansV = 0
    for (nombreDeDimension in nbreVarX) {
        nbreVarX_dansV = nbreVarX_dansV + 1
        for (b in 1:length(Y)) {
            apprends = matrix(X[-b, ], ncol = dim(X)[2])
            verif = matrix(X[b, ], ncol = dim(X)[2])
            param.fs[["X"]] = apprends
            param.fs[["Y"]] = Y[-b]
            param.fs[["nbreVarX"]] = nombreDeDimension
            f = do.call(filter, param.fs)
            param.model[["x"]] = matrix(predict(f, apprends)[,1: nombreDeDimension],ncol=nombreDeDimension)
            param.model[["y"]] = factor(Y[-b])
            predicteur = do.call(model, param.model)
            p <- predict(predicteur, c(predict(f, verif)))
            if (length(unique(Y))==2) {
                if (as.numeric(levels(p))[p]==0 && Y[b]==0) {
                    TN[nbreVarX_dansV] <- TN[nbreVarX_dansV] + 1
                }
                if (as.numeric(levels(p))[p]==0 && Y[b]!=0) {
                    FN[nbreVarX_dansV] <- FN[nbreVarX_dansV] + 1
                }
                if (as.numeric(levels(p))[p]!=0 && Y[b]==0) {
                    FP[nbreVarX_dansV] <- FP[nbreVarX_dansV] + 1
                }
                if (as.numeric(levels(p))[p]!=0 && Y[b]!=0) {
                    TP[nbreVarX_dansV] <- TP[nbreVarX_dansV] + 1
                }
            }
            else{
                if(as.numeric(levels(p))[p]==Y[b]) {
                    TN[nbreVarX_dansV] <- TN[nbreVarX_dansV] + 1
                }
                else {
                    FN[nbreVarX_dansV] <- FN[nbreVarX_dansV] + 1
                }
            }
        }
    }
    res = list(TP = TP, TN = TN, FN = FN, FP = FP, dim = nbreVarX)
    class(res) <- "verify.loo"
    return(res)
}
verify.ho.default <-
function (model, filter, Xlearn, Ylearn, Xval, Yval, nbreVarX, param.model=list(), param.fs=list(), ...) 
{
    if (!is.vector(Ylearn) || !is.vector(Yval)) {
        stop("'Y' has to be a vector")
        return(-1)
    }
    if (!is.matrix(Xlearn) || !is.matrix(Xval)) {
        stop("'X' has to be a matrix")
        return(-1)
    }
    if (( length(Ylearn) != dim(Xlearn)[1]) || ( length(Yval) != dim(Xval)[1])) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
    TP <- rep(0, length(nbreVarX))
    TN <- rep(0, length(nbreVarX))
    FN <- rep(0, length(nbreVarX))
    FP <- rep(0, length(nbreVarX))
    nbreVarX_dansV = 0
    for (nombreDeDimension in nbreVarX) {
        nbreVarX_dansV = nbreVarX_dansV + 1
        param.fs[["X"]] = Xlearn
        param.fs[["Y"]] = Ylearn
        param.fs[["nbreVarX"]] = nombreDeDimension
        f = do.call(filter, param.fs)
        param.model[["x"]] = matrix(predict(f, Xlearn),ncol=nombreDeDimension)
        param.model[["y"]] = factor(Ylearn)
        predicteur = do.call(model, param.model)
	p <- predict(predicteur, matrix(predict(f, Xval),ncol=nombreDeDimension))
	p <- as.numeric(levels(p))[p]
	if (length(unique(Ylearn))==2) {
		TN[nbreVarX_dansV] <- TN[nbreVarX_dansV] + length( which(p==0 & Yval==0))
        	FN[nbreVarX_dansV] <- FN[nbreVarX_dansV] + length( which(p==0 & Yval==1))
        	FP[nbreVarX_dansV] <- FP[nbreVarX_dansV] + length( which(p==1 & Yval==0))
        	TP[nbreVarX_dansV] <- TP[nbreVarX_dansV] + length( which(p==1 & Yval==1))
	}
	else{
		TN[nbreVarX_dansV] <- TN[nbreVarX_dansV] + length( which(p==Yval))
        	FN[nbreVarX_dansV] <- FN[nbreVarX_dansV] + length( which(p!=Yval))
	}
    }
    res = list(TP = TP, TN = TN, FN = FN, FP = FP, dim = nbreVarX)
    class(res) <- "verify.ho"
    return(res)
}
verify.cv.default <-
function (model, filter, X, Y, nbreVarX, k, param.model=list(), param.fs=list(), ...) 
{
    if (!is.vector(Y)) {
        stop("'Y' has to be a vector")
        return(-1)
    }
    if (!is.matrix(X)) {
        stop("'X' has to be a matrix")
        return(-1)
    }
    if (length(Y) != dim(X)[1]) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
    if (k >= dim(X)[1]) {
        stop("the number k has to be less than the number of input")
        return(-1)
    }


    TP <- rep(0, length(nbreVarX))
    TN <- rep(0, length(nbreVarX))
    FN <- rep(0, length(nbreVarX))
    FP <- rep(0, length(nbreVarX))
    divK = seq(1,dim(X)[1],k)
    for(b in 1:(length(divK)-1)){
        apprendsX = matrix(X[-c(divK[b]:(divK[b+1]-1)), ], ncol = dim(X)[2])
        verifX = matrix(X[c(divK[b]:(divK[b+1]-1)), ], ncol = dim(X)[2])
        apprendsY = Y[-c(divK[b]:(divK[b+1]-1)) ]
        verifY = Y[c(divK[b]:(divK[b+1]-1)) ]
        res = verify.ho(model=model, filter=filter, Xlearn=apprendsX, Ylearn=apprendsY, Xval=verifX, Yval=verifY, nbreVarX=nbreVarX, param.model, param.fs)
        TP = TP + res$TP
        TN = TN + res$TN
        FN = FN + res$FN
        FP = FP + res$FP
    }
    apprendsX = matrix(X[-c(divK[length(divK)]:(dim(X)[1])), ], ncol = dim(X)[2])
    verifX = matrix(X[c(divK[length(divK)]:(dim(X)[1])), ], ncol = dim(X)[2])
    apprendsY = Y[-c(divK[length(divK)]:(dim(X)[1])) ]
    verifY = Y[c(divK[length(divK)]:(dim(X)[1])) ]
    res = verify.ho(model=model, filter=filter, Xlearn=apprendsX, Ylearn=apprendsY, Xval=verifX, Yval=verifY, nbreVarX=nbreVarX, param.model, param.fs)
    TP = TP + res$TP
    TN = TN + res$TN
    FN = FN + res$FN
    FP = FP + res$FP
    res = list(TP = TP, TN = TN, FN = FN, FP = FP, dim = nbreVarX)
    class(res) <- "verify.loo"
    return(res)
}
