benchmark.regression=function (X, y, m = ncol(X), R = 20, ratio = 0.8, verbose = TRUE,k = 10, nsamples = nrow(X), use.kernel = FALSE,supervised=FALSE) {
    n <- nsamples
    m.pls <- m.pcr<-lambda.ridge<-vector(length = R) # vector of optimal model parameters
    ntrain <- floor(n * ratio) # number of training observations
    DoF.pls <- DoF.pcr <- DoF.ridge<-vector(length = R) # Degrees of Freedom
    mse.pls <- mse.pcr <- mse.ridge <- vector(length = R) # test error
    mse.null <- vector(length = R) # test error of the null model, i.e. mean(ytrain)
    DoF.all<-res.pls<-res.pcr<-matrix(,R,ncol(X)+1) # residuals, do we need this, do we need DoF.all
    for (i in 1:R) {
        if (verbose == TRUE) {
            cat(paste("iteration no ", i, " \n"))
        }
        # subsample
        samples <- sample(nrow(X), n, replace = FALSE)
        XX <- X[samples, ]
        yy <- y[samples]
        # split into training and test
        train <- sample(n, ntrain, replace = FALSE)
        Xtrain <- XX[train, , drop = FALSE]
        Xtest <- XX[-train, , drop = FALSE]
        ytrain <- yy[train]
        ytest <- yy[-train]
        # null model
        mse.null[i] = mean((mean(ytrain) - ytest)^2)
        # pls
        pls.object<- pls.cv(Xtrain, ytrain, use.kernel = use.kernel, m = m,k=k)
        m.pls[i]<- pls.object$m.opt
        pls.object = pls.model(Xtrain, ytrain, Xtest = Xtest,ytest = ytest, m=m,compute.DoF = TRUE, compute.jacobian = FALSE, use.kernel = use.kernel)
        res.train<-ytrain-rep(1,ntrain)%*%t(pls.object$intercept) - Xtrain%*%pls.object$coefficients
        res.pls[i,]<-apply(res.train^2,2,mean)
        DoF.all[i,]<-pls.object$DoF
        mse.pls[i] <- pls.object$mse[m.pls[i] + 1]
        DoF.pls[i]<-pls.object$DoF[m.pls[i]+1]
        #pcr
        pcr.object<-pcr.cv(Xtrain,ytrain,k=k,m=m,supervised=supervised)
        m.pcr[i]<-pcr.object$m.opt
        DoF.pcr[i]<-m.pcr[i]+1
        pcr.final<-pcr(Xtrain,ytrain)
        res.train<-ytrain - rep(1,ntrain)%*%t(pcr.final$intercept) - Xtrain%*%pcr.final$coefficients
        res.pcr[i,]<-apply(res.train^2,2,mean)
        res<-ytest-pcr.object$intercept - Xtest%*%pcr.object$coefficients
        mse.pcr[i]<-apply(res^2,2,mean)
        # ridge
        ridge.object<-ridge.cv(Xtrain,ytrain,k=k)
        lambda.ridge[i]<-ridge.object$lambda.opt
        res<-ytest-ridge.object$intercept - Xtest%*%ridge.object$coefficients
        mse.ridge[i]<-apply(res^2,2,mean)
        XX<-scale(Xtrain)
        SS<-t(XX)%*%XX
        lambda<-eigen(SS)$values
        DoF.ridge[i]=sum(lambda/(lambda+lambda.ridge[i]))
        }
    namen <- c("PLS", "PCR", "RIDGE","NULL")
    MSE <- matrix(, R, 4)
    MSE[, 1] <- mse.pls
    MSE[, 2] <- mse.pcr
    MSE[, 3] <- mse.ridge
    MSE[,4]<-mse.null
    MSE <- data.frame(MSE)
    colnames(MSE) = namen
    M <- matrix(, R, 4)
    M[, 1] <- m.pls
    M[, 2] <- m.pcr
    M[, 3] <- lambda.ridge
    M[, 4] <- rep(0, R)
    M <- data.frame(M)
    colnames(M) = namen
    DoF <- matrix(, R, 4)
    DoF[, 1] <- DoF.pls
    DoF[, 2] <- DoF.pcr
    DoF[, 3] <- DoF.ridge
    DoF[, 4] <- rep(1, R)
    DoF <- data.frame(DoF)
    colnames(DoF) = namen
    return(list(MSE = MSE, M = M, DoF = DoF,res.pls=res.pls,res.pcr=res.pcr,DoF.all=DoF.all))
}
