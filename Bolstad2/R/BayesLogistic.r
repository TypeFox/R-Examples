BayesLogistic = function(y, x, steps = 1000,
                        priorMean = NULL, priorVar = NULL,
                        mleMean = NULL, mleVar,
                        startValue = NULL, randomSeed = NULL,
                        plots = FALSE){

    if(!is.null(randomSeed))
        set.seed(randomSeed)

    nObs = length(y)

    if(is.vector(x))
        x = as.matrix(x, ncol = 1)

    nParameters = ncol(x) + 1 ## number of covariates + intercept

    if(!is.null(startValue)){
        if(length(startValue) < nParameters){
            stop("You must have as many starting values as you have model parameters")
        }
    }

    ## inital mean of the matched curvature likelihood
    if(is.null(mleMean))
        mleMean = c(log(mean(y)/(1-mean(y))), rep(0, nParameters-1))

    X = cbind(rep(1,nObs),x)
    Xt = t(X)


    calcMatchedCurvatureNormLike = function(){

        betaX = X%*%mleMean
        Pi = exp(betaX)/(1+exp(betaX))
        Vdiag = Pi*(1-Pi)
        Y = betaX + (y-Pi)/Vdiag

        ## I have no idea why the diag command doesn't work as it should:
        ## e.g.            Vyinv = diag(Vdiag, nrow = length(Vdiag))
        ## therefore this two-step procedure is needed
        Vyinv = matrix(0, nrow = nObs, ncol = nObs)
        diag(Vyinv) = Vdiag

        XtV = Xt%*%Vyinv
        VLinv = XtV%*%X
        VL = solve(VLinv)
        w1 = VL%*%XtV
        mleMean = w1%*%Y


        ##   Loop iterations to converge to MLE
        for(k in 1:5){
            betaX = X%*%mleMean
            Pi = exp(betaX)/(1+exp(betaX))
            Vdiag = Pi*(1-Pi)
            Y = betaX + (y-Pi)/Vdiag
            Vyinv = matrix(0, nrow = nObs, ncol = nObs)
            diag(Vyinv) = Vdiag

            XtV = Xt%*%Vyinv
            VLinv = XtV%*%X
            VL = solve(VLinv)
            w1 = VL%*%XtV
            mleMean = w1%*%Y
        }

        return(list(mleMean = mleMean, mleVar = VL))
    } ## calcMatchedCurvatureNormLike

    normApproxPosterior = function(){

        result = list(postMean = rep(0, nParameters),
                       postVar = matrix(0, ncol = nParameters, nrow = nParameters))



        ## if the prior mean and variance isn't specified then
        ## set it equal to the mle mean and variance
        if(is.null(priorMean)){
            priorMean =  result$postMean = mleMean
        }

        if(is.null(priorVar)){
            priorVar = result$postVar = mleVar
        }

        mleVarInv = solve(mleVar)
        priorVarInv = solve(priorVar)
        postPrec = mleVarInv + priorVarInv
        result$postVar = solve(postPrec)

        result$postMean = result$postVar%*%priorVarInv%*%priorMean + result$postVar%*%mleVarInv%*%mleMean

        return(result)
    }

    mleParams = calcMatchedCurvatureNormLike()
    mleMean = mleParams$mleMean
    mleVar = mleParams$mleVar

    posterior = normApproxPosterior()
    postMean = posterior$postMean
    postVar = posterior$postVar

    U = chol(postVar)
    L = t(U)

    candBeta = matrix(rt(steps*nParameters, df = 4), ncol = nParameters)

    if(!is.null(startValue))
        candBeta[1,] = startValue

    WM2 = candBeta %*% U
    WM3 = matrix(rep(postMean,rep(steps,nParameters)),ncol = nParameters)
    WM4 = WM2 + WM3
    V2 = cov(WM4)

    ft0 = apply(dt(candBeta, df = 4), 1, prod)
    fn0 = apply(dnorm(candBeta), 1, prod)
    q1 = ft0/1

    ## Metropolis-Hastings

    Sum1 = WM4 %*% Xt

    Pi1 = exp(Sum1)/(1+exp(Sum1))

    for(j in 1:nObs)
        Pi1[,j] = log(Pi1[,j]^y[j]*(1-Pi1[,j])^(1-y[j]))

    g0 = exp(rowSums(Pi1))
    g1 = g0

    if(!is.null(priorMean))
        g1 = g0*fn0

    g1 = g1/max(g1)
    q1 = q1/max(q1)

    plot(q1, g1)

    u = runif(steps)
    i1 = 1

    betaSample = WM4

    for(n in 2:steps){
        alpha = q1[i1]*g1[n]/(q1[n]*g1[i1])
        alpha = ifelse(alpha>1, 1, alpha)

        if(u[n] >= alpha){ ## reject
            betaSample[n,] = WM4[i1,]
        }else{
            betaSample[n,] = WM4[n,]
            i1 = n
        }
    }

    beta.df = data.frame(betaSample)
    names(beta.df) = paste("b",0:(ncol(beta.df)-1),sep = "")
    describe(beta.df)

    Mean.beta = sapply(beta.df,mean)
    StdDev.beta = sapply(beta.df,sd)
    Z.beta = Mean.beta/StdDev.beta

    print(data.frame(Mean.beta,StdDev.beta,Z.beta))

    if(plots){
     ##   nRows = ceiling(sqrt(nParameters))
        nRows = nParameters
     ##   nCols = floor(sqrt(nParamerts))
        nCols = 2
        oldPar = par(mfrow = c(nRows, nCols))
        nms = names(beta.df)

        for(i in 1:nParameters){
            plot(ts(beta.df[,i]),
                 main = paste("Time series plot of",nms[i]),
                 ylab = nms[i])
            plot(acf(beta.df[,i], plot = FALSE),
                 main = paste("Autocorrelation plot of", nms[i]))
        }

        par(oldPar)
    }

    invisible(list(beta = beta.df, mleMean = mleMean, mleVar = mleVar))
}
