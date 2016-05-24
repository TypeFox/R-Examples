auc.para.bayes <- function(x, y, conf.level=0.95, 
                           dist=c("normalDV", "normalEV", "exponential"),
                           nburn=1000, niter=10000, nthin=1, seed=100){
    
    dist <- match.arg(dist)
    R <- switch(dist,
    normalEV={
        #two normal with different means and same variance
        modelstring = "
        model {
        for(i in 1:nx){
            x[i] ~ dnorm (mux, tau) 
        }
        for(i in 1:ny){
            y[i] ~ dnorm (muy, tau) 
        }
        mux ~ dnorm (0, 0.0001) 
        muy ~ dnorm (0, 0.0001) 
        tau ~ dgamma (0.001, 0.001) 
        }
        "
        writeLines(modelstring, con="model.txt")

        nx <- length(x)
        ny <- length(y)
        dataList <- list(x=x, y=y, nx=nx, ny=ny)
        x.s2 <- var(x)
        y.s2 <- var(y)
        s2 <- ((nx-1)*x.s2 + (ny-1)*y.s2) / (nx+ny-2) 
        initsList <- list(mux=mean(x), muy=mean(y), tau=1/s2, 
                          .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
        parameters = c("mux", "muy", "tau")
        jagsModel = jags.model("model.txt" , data=dataList , inits=initsList , 
                               n.chains=1, n.adapt=0, quiet=TRUE)
        
        update(jagsModel , n.iter=nburn, progress.bar="none")
        codaSamples = coda.samples(jagsModel, variable.names=parameters, 
                                   n.iter=niter, thin=nthin, progress.bar="none")
        
        mcmcChain = as.matrix(codaSamples)
        mux <- mcmcChain[,1]
        muy <- mcmcChain[,2]
        sigma2 <- 1/mcmcChain[,3]
        pnorm((mux- muy)/sqrt(2*sigma2))
        
    },
    normalDV={
        modelstring = "
        #two normal with different means and different variances
        model {
        for(i in 1:nx){
            x[i] ~ dnorm (mux, taux) 
        }
        for(i in 1:ny){
            y[i] ~ dnorm (muy, tauy) 
        }
        mux ~ dnorm (0, 0.0001) 
        muy ~ dnorm (0, 0.0001) 
        taux ~ dgamma (0.001, 0.001) 
        tauy ~ dgamma (0.001, 0.001) 
        }
        " 
        writeLines(modelstring, con="model.txt")
        
        dataList <- list(x=x, y=y, nx=length(x), ny=length(y))
        x.s2 <- var(x)
        y.s2 <- var(y)
        initsList <- list(mux=mean(x), muy=mean(y), taux=1/x.s2, tauy=1/y.s2,
                          .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
        parameters <- c("mux", "muy", "taux", "tauy")
        jagsModel = jags.model("model.txt" , data=dataList , inits=initsList , 
                               n.chains=1, n.adapt=0, quiet=TRUE)
        
        update(jagsModel , n.iter=nburn, progress.bar="none")
        codaSamples = coda.samples(jagsModel, variable.names=parameters, 
                                   n.iter=niter, thin=nthin, progress.bar="none")
        
        mcmcChain  <-  as.matrix(codaSamples)
        mux <- mcmcChain[,1]
        muy <- mcmcChain[,2]
        sigma2x <- 1 / mcmcChain[,3]
        sigma2y <- 1 / mcmcChain[,4]
        pnorm((mux- muy)/sqrt(sigma2x+sigma2y))
        
    },
    exponential={
        modelstring = "
        #two normal with different means and different variances
        model {
        for(i in 1:nx){
            x[i] ~ dexp (lambdax) 
        }
        for(i in 1:ny){
            y[i] ~ dexp (lambday) 
        }
        lambdax ~ dgamma (0.001, 0.001) 
        lambday ~ dgamma (0.001, 0.001) 
        }
        " 
        writeLines(modelstring, con="model.txt")
        
        dataList <- list(x=x, y=y, nx=length(x), ny=length(y))
        initsList <- list(lambdax=1/mean(x), lambday=1/mean(y),
                          .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
        parameters <- c("lambdax", "lambday")
        jagsModel = jags.model("model.txt" , data=dataList , inits=initsList , 
                               n.chains=1, n.adapt=0, quiet=TRUE)
        
        update(jagsModel , n.iter=nburn, progress.bar="none")
        codaSamples = coda.samples(jagsModel, variable.names=parameters, 
                                   n.iter=niter, thin=nthin, progress.bar="none")
        
        mcmcChain  <-  as.matrix(codaSamples)
        lambdax <- mcmcChain[,1]
        lambday <- mcmcChain[,2]
        lambday / (lambdax + lambday)
        
    })
    
    c(mean(R), as.vector(HPDinterval(mcmc(R, conf.level))))
}
