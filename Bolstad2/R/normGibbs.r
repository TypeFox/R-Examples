normGibbs = function(y, steps = 1000, type = 'ind', ...){

    if(length(grep('[Ii]',type))>0){
        type = 'ind'
    }else if(length(grep('[Jj]',type))>0){
        type = 'joint'
    }else{
        stop("Type must be ind or joint")
    }

    dots = list(...)


    if(type == 'ind'){

        ## The dots can carry priorMu which, if specified
        ## can be a single number m0 or a vector m0, n0

        m0 = 0
        n0 = 0
        muSpecified = FALSE

        if("priorMu" %in% names(dots)){
            if(length(dots$priorMu)==2){
                m0 = dots$priorMu[1]
                n0 = dots$priorMu[2]
                muSpecified = TRUE
            }else{
                m0 = dots$priorMu
                muSpecified = TRUE
            }
        }

        ## The dots can carry priorVar which, if specified
        ## must be a vector of length 2

        s0 = 0
        kappa0 = 0
        nObs = length(y)
        yBar = mean(y)
        SSy = sum((y-yBar)^2)

        varSpecified = FALSE

        if("priorVar" %in% names(dots)){
            if(length(dots$priorVar)!=2)
                stop("priorVar must have two elements, s0 and kappa0")
            else{
                s0 = dots$priorVar[1]
                kappa0 =dots$priorVar[2]
                varSpecified = TRUE
            }
        }

        v0 = ifelse(varSpecified, s0/rchisq(1,df = kappa0), SSy/(nObs-1))
        kappa1 = kappa0 + nObs

        prec0 = 0
        mu0 = yBar

        if(muSpecified){
            prec0 = n0/v0
            mu0 = m0 + rnorm(1)*sqrt(v0)
        }else{
            prec0 = 0
            m0 = 0
            mu0 = yBar
        }

        SSt = sum((y-mu0)^2)
        s1 = s0 + SSt

        Chi = rchisq(steps, df = kappa1)
        z = rnorm(steps)

        varSample = c(s1/Chi[1],rep(0,steps-1))


        prec0 = ifelse(muSpecified, n0/varSample[1], 0)

        precData = nObs / varSample[1]
        prec1 = prec0 + precData
        v1 = 1/prec1
        m1 = m0*prec0/prec1 + yBar*precData/prec1

        muSample = rep(0, steps)
        muSample[1] = z[1]*sqrt(v1)+m1

        for(i in 2:steps){
            SSt = sum((y-muSample[i-1])^2)
            s1 = s0 + SSt
            varSample[i] = s1/Chi[i]

            prec0 = ifelse(muSpecified, n0/varSample[i], 0)
            precData = nObs / varSample[i]
            prec1 = prec0 + precData
            v1 = 1/prec1
            std1 = sqrt(v1)
            m1 = m0*prec0/prec1 + yBar*precData/prec1

            muSample[i] = z[i]*std1 + m1
        }

        sigmaSample = sqrt(varSample)

        oldPar = par(mfrow = c(2,2))

        plot(ts(muSample),main = "Time series plot of mu", ylab = expression(mu))
        plot(ts(varSample),main = "Time series plot of var"
             , ylab = expression(sigma^2))
        plot(ts(sigmaSample), main = "Time series plot of sigma",
             ylab = expression(sigma))
        hist(muSample, main = "Histogram of mu", prob = TRUE)
        mx = mean(muSample)
        sx = sd(muSample)
        bds = mx+c(-5,5)*sx
        xValues=seq(bds[1],bds[2],length = 200)
        yValues = dnorm(xValues,mx,sx)
        lines(xValues,yValues)

        par(oldPar)

        results.df = data.frame(mu = muSample, sig = sigmaSample, var = varSample)

        describe(results.df)
        invisible(results.df)
    }else{ ## type = 'joint'

        ## The dots can carry priorMu which, if specified
        ## can be a single number m0 or a vector m0, n0

        m0 = 0
        n0 = 0
        v0 = 0

        muSpecified = FALSE

        if("priorMu" %in% names(dots)){
            if(length(dots$priorMu)==2){
                m0 = dots$priorMu[1]
                n0 = dots$priorMu[2]
                muSpecified = TRUE
            }else{
                stop("priorMu must contain a mean and an effective sample size")
            }
        }

        ## The dots can carry priorVar which, if specified
        ## must be a vector of length 2

        s0 = 0
        kappa0 = 0
        varSpecified = FALSE

        v0 = 0

        if("priorVar" %in% names(dots)){
            if(length(dots$priorVar)!=2)
                stop("priorVar must have two elements, s0 and kappa0")
            else{
                s0 = dots$priorVar[1]
                kappa0 =dots$priorVar[2]
                varSpecified = TRUE
            }
        }

        nObs = length(y)
        yBar = mean(y)
        SSy = sum((y-yBar)^2)

        v0 = ifelse(varSpecified, s0/rchisq(1,df = kappa0), SSy/(nObs-1))

        kappa1 = kappa0 + nObs

        prec0 = 0
        mu0 = yBar

        muSample = rep(0, steps)

        if(muSpecified){
            prec0 = n0/v0
            mu0 = m0 + rnorm(1)*sqrt(v0)
        }else{
            m0 = 0
        }


        SSt = sum((y-mu0)^2)
        s1 = s0 + SSt

        Chi = rchisq(steps, df = kappa1)

        varSample = c(s1/Chi[1],rep(0,steps-1))
        Z = rnorm(steps)

        prec0 = ifelse(muSpecified, n0/varSample[1], 0)
        precData = nObs/varSample[1]
        prec1 = prec0 + precData
        v1 = 1/prec1
        m1 = prec0/prec1*m0 + precData/prec1*yBar

        muSample[1] = Z[1]*sqrt(v1)+m1

        for(i in 2:steps){
            SSt = sum((y-muSample[i-1])^2)
            s1 = s0 + SSt
            varSample[i] = s1/Chi[i]

            prec0 = ifelse(muSpecified, n0/varSample[i],0)
            precData = nObs/varSample[i]
            prec1 = prec0 + precData
            v1 = 1/prec1
            m1 = prec0/prec1*m0 + precData/prec1*yBar

            muSample[i] = Z[i]*sqrt(v1)+m1
        }

        sigmaSample = sqrt(varSample)

        oldPar = par(mfrow = c(2,2))

        plot(ts(muSample),main = "Time series plot of mu", ylab = expression(mu))
        plot(ts(varSample),main = "Time series plot of var"
             , ylab = expression(sigma^2))
        plot(ts(sigmaSample), main = "Time series plot of sigma",
             ylab = expression(sigma))
        hist(muSample, main = "Histogram of mu", prob = TRUE)
        mx = mean(muSample)
        sx = sd(muSample)
        bds = mx+c(-5,5)*sx
        xValues=seq(bds[1],bds[2],length = 200)
        yValues = dnorm(xValues,mx,sx)
        lines(xValues,yValues)

        par(oldPar)

        results.df = data.frame(mu = muSample, sig = sigmaSample, var = varSample)

        describe(results.df)
        invisible(results.df)
    }
}





