##' simsurv function
##'
##' A function to simulate spatial parametric proportional hazards model. The function works
##' by simulating candidate survival times using MCMC in parallel for each individual based on each individual's covariates and the common
##' parameter effects, beta.  
##'
##' @param X a matrix of covariate information 
##' @param beta the parameter effects 
##' @param omega vector of parameters for the baseline hazard model
##' @param dist the distribution choice: exp or weibull at present
##' @param coords matrix with 2 columns giving the coordinates at which to simulate data
##' @param cov.parameters a vector: the parameters for the covariance function
##' @param cov.model an object of class covmodel, see ?covmodel
##' @param mcmc.control mcmc control paramters, see ?mcmcpars
##' @param savechains save all chains? runs faster if set to FALSE, but then you'll be unable to conduct convergence/mixing diagnostics
##' @return in list element 'survtimes', a vector of simulated survival times (the last simulated value from the MCMC chains) 
##' in list element 'T' the MCMC chains
##' @seealso \link{covmodel}, \link{survspat}, \link{tpowHaz}, \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz} 
##' @export

simsurv <- function(X=cbind(age=runif(100,5,50),sex=rbinom(100,1,0.5),cancer=rbinom(100,1,0.2)),
                            beta=c(0.0296,0.0261,0.035),
                            omega=1,
                            dist="exp",
                            coords=matrix(runif(2*nrow(X)),nrow(X),2),
                            cov.parameters=c(1,0.1),
                            cov.model=covmodel(model="exponential",pars=NULL),
                            mcmc.control=mcmcpars(nits=100000,burn=10000,thin=90),      
                            savechains=TRUE){                         

    beta <- matrix(beta,length(beta),1)

    mcmcloop <- mcmcLoop(N=mcmc.control$nits,burnin=mcmc.control$burn,thin=mcmc.control$thin,progressor=mcmcProgressTextBar)

    distmat <- as.matrix(stats::dist(coords))

    n <- nrow(X)
    u <- as.vector(distmat)

    sigma <- matrix(EvalCov(cov.model,u=u,parameters=cov.parameters),n,n)  
    
    sigmachol <- t(chol(sigma))
    Y <- -cov.parameters[which(cov.model$parnames=="sigma")]^2/2 + sigmachol%*%rnorm(n)
    expY <- exp(Y)  
 
    XbetaplusY <- X%*%beta + Y
    expXbetaplusY <- exp(XbetaplusY)
    
    h <- basehazard(dist)(omega)    
    H <- cumbasehazard(dist)(omega)
    
    nmatrows <- ceiling((mcmc.control$nits-mcmc.control$burn-mcmcloop$waste)/mcmc.control$thin)
    
    if(savechains){
        T <- matrix(NA,nmatrows,n) #matrix to store simulated times
    }
    else{
        T <- NULL
    }
    tarrec <- matrix(NA,nmatrows,n)
    acrec <- rep(NA,nmatrows)
    count <- 1 # counter to index retained iteration numbers
    
    #if(dist=="exp"){
    #    t <- rexp(n,omega)
    #}
    #else if(dist=="weibull"){
    #    transpars <- transformweibull(omega)
    #    t <- rweibull(n,shape=transpars[1],scale=transpars[2])
    #}
    
    t <- rep(1,n)
    
    oldltar <- XbetaplusY + log(h(t)) - expXbetaplusY*H(t) 
    
    tbar <- rep(0,n)
    nsamp <- 0
    
    while(nextStep(mcmcloop)){        
    
        newt <- rexp(n,rate=1/t) 
        newltar <- XbetaplusY + log(h(newt)) - expXbetaplusY*H(newt)
        
        frac <- exp(newltar-oldltar+dexp(t,1/newt,log=TRUE)-dexp(newt,1/t,log=TRUE))        
        
        ac <- pmin(1,frac)
        
        keepnew <- ac>runif(n)
        t[keepnew] <- newt[keepnew]
        oldltar[keepnew] <- newltar[keepnew]
    
        if(is.retain(mcmcloop)){
            nsamp <- nsamp + 1 # note the purpose of this counter is different to the "count" counter
            if(savechains){        
                T[count,] <- t
            }
            
            tbar <- ((nsamp-1)/nsamp)*tbar + (1/nsamp)*t            
            
            tarrec[count,] <- oldltar
            acrec[count] <- mean(ac)
            count <- count + 1
        }        
    }   
    
    if(is.na(tarrec[nmatrows])){
        if(savechains){
            T <- T[-nmatrows,]
        }
        tarrec <- tarrec[-nmatrows]
        acrec <- acrec[-nmatrows]
    }
    
    cat("Returning last set of simulated survival times.\n")
    tchoice <- t
    
    cat("Mean acceptance:",mean(acrec),"\n")
    
    return(list(X=X,T=T,survtimes=tchoice,omega=omega,beta=beta,tarrec=tarrec,n=n,Y=Y,dist=dist,coords=coords,cov.parameters=cov.parameters,distmat=distmat,u=u))
}