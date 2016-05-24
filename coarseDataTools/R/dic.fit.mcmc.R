##' Fits the distribution to the passed-in data using MCMC
##' as implemented in MCMCpack.
##'
##' Similar to \code{\link{dic.fit}} but uses MCMC instead of a direct likelihood optimization routine to fit the model. Currently, four distributions are supported: log-normal, gamma, Weibull, and Erlang. See Details for prior specification.
##'
##' The following models are used:
##' \deqn{Log-normal model: f(x) = \frac{1}{x*\sigma \sqrt{2 * \pi}} exp\{-\frac{(\log x - \mu)^2}{2 * \sigma^2}\}}
##' \deqn{Log-normal Default Prior: \mu ~ N(0, 1000), log(\sigma) ~ N(0,1000)}
##' \deqn{Weibull model: f(x) = \frac{\alpha}{\beta}(\frac{x}{\beta})^{\alpha-1} exp\{-(\frac{x}{\beta})^{\alpha}\}}
##' \deqn{Weibull Default Prior Specification: log(\alpha) ~ N( 0, 1000), \beta ~ Gamma(0.001,0.001)}
##' \deqn{Gamma model: f(x) = \frac{1}{\theta^k \Gamma(k)} x^{k-1} exp\{-\frac{x}{\theta}\}}
##' 
##' \deqn{Gamma Default Prior Specification: p(k,\theta) \propto \frac{1}{\theta} * \sqrt{k*TriGamma(k)-1}}
##' (Note: this is Jeffery's Prior when both parameters are unknown), and 
##' \deqn{Trigamma(x) = \frac{\partial}{\partial x^2} ln(\Gamma(x))}.)
##' \deqn{Erlang model: f(x) = \frac{1}{\theta^k (k-1)!} x^{k-1} exp\{-\frac{x}{\theta}\}}
##' \deqn{Erlang Default Prior Specification: k \sim NBinom(100,1), log(\theta) \sim N(0,1000)}
##' (Note: parameters in the negative binomial distribution above represent mean and size, respectivley)
##'
##' @param dat the data
##' @param prior.par1 vector of first prior parameters for each model parameter. If \code{NULL} then default parameters are used (as described in Details section).
##' @param prior.par2 vector of second prior parameters for each model parameter. If \code{NULL} then default parameters are used (as described in Details section).
##' @param init.pars the initial parameter values (vector length = 2 )
##' @param ptiles returned percentiles of the survival survival distrbution
##' @param verbose how often do you want a print out from MCMCpack on iteration number and M-H acceptance rate
##' @param burnin number of burnin samples
##' @param n.samples number of samples to draw from the posterior (after the burnin)
##' @param dist distribution to be used (L for log-normal,W for weibull, G for Gamma, and E for erlang)
##' @param ... additional parameters to \link{MCMCmetrop1R}
##' @return a cd.fit.mcmc S4 object
##' @importFrom MCMCpack MCMCmetrop1R
##' @export
dic.fit.mcmc <- function(dat,
                         prior.par1 = NULL,
                         prior.par2 = NULL,
                         init.pars = c(1,1),
                         ptiles = c(0.05,0.95,0.99),
                         verbose=1000,#how often to print update
                         burnin = 3000,
                         n.samples = 5000,
                         dist = "L",
                         ...){

        
        ## check to make sure data is well formed for CDT use:
        check.data.structure(dat)
                       
        ## check to make sure distribution is supported
        if(!dist %in% c("G","W","L","E")) stop("Please use one of the following distributions Log-Normal (L) , Weibull (W), Gamma (G), or Erlang (E)")
        
##        ## don't need MCMCpack for Erlang
##        if (dist != "E") require(MCMCpack)       
         
        ## check that percentiles are valid
        if (any(ptiles >=1) | any(ptiles <= 0)) stop("Sorry the percentiles you are requesting are not valid.")
        
        ## default prior parameters if none specified
        if (is.null(prior.par1)){
           if (dist == "L") {
            prior.par1 <- c(0,0)
           } else if (dist == "W" | dist == "G"){
            prior.par1 <- c(0,0.001)
           } else if (dist == "E"){
            prior.par1 <- c(100,0)          
         }}
        
        ## default prior parameters if none specified        
        if (is.null(prior.par2)){
         if (dist == "L") {
          prior.par2 <- c(1000,1000)
         } else if (dist == "W" | dist == "G"){
          prior.par2 <- c(1000,0.001)
         } else if (dist == "E"){
          prior.par2 <- c(1,1000)          
         }
        }
              
        cat(sprintf("Running %.0f MCMC iterations \n",n.samples+burnin))

        msg <- NULL
        fail <- FALSE

        ## run the MCMC chains

        ## initial parameters set on reporting scale not estimation scale
        init.pars.trans <- dist.optim.transform(dist,init.pars)
  
      
        if (dist!="E") {
          #use MCMC pack for effieciency for distibutions with 2 continuous parameters
          tryCatch(            
              mcmc.run <- MCMCmetrop1R(fun=mcmcpack.ll,
                                          theta.init=init.pars.trans,
                                          burnin = burnin,
                                          mcmc=n.samples,
                                          dat = dat,
                                          prior.par1 = prior.par1,
                                          prior.par2 = prior.par2,
                                          verbose=verbose,
                                          dist=dist,
                                          logfun=TRUE,
                                          ...),
                 error=function(e){
                         msg <<- e$message
                         fail <<- TRUE
                 },
                 warning = function(w){
                         msg <<- w$message
                         fail <<- TRUE
                 })
        } else {
          #Use the custom metropolis hastings algorithm for erlang
          mcmc.run <- mcmc.erlang(dat, 
                                  prior.par1,
                                  prior.par2,
                                  init.pars, 
                                  verbose,
                                  burnin,
                                  n.samples,
                                  ...)
        }
        
        

        if (!fail){
                ## untransform MCMC parameter draws to natural scale
                untrans.mcmcs <- t(apply(mcmc.run[,1:2],1,function(x) dist.optim.untransform(dist=dist,pars=x)))

                ## append median to the percentiles in case it isn't there
                ptiles.appended <- sort(union(0.5,ptiles))
                est.pars <- matrix(nrow=length(ptiles.appended)+2,ncol=3)

                if (dist == "L"){
                        par1.name <- "meanlog"
                        par2.name <- "sdlog"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qlnorm(ptiles.appended,meanlog=x[1],sdlog=x[2]))
                } else if (dist == "G"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))
                } else if (dist == "W"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qweibull(ptiles.appended,shape=x[1],scale=x[2]))
                } else if (dist == "E"){
                        par1.name <- "shape"
                        par2.name <- "scale"
                        mcmc.quantiles <- apply(untrans.mcmcs,1,function(x) qgamma(ptiles.appended,shape=x[1],scale=x[2]))
                } else {
                        stop("Sorry, unknown distribution type. Check the 'dist' option.")
                        ## not actually needed but just in case
                }

                ## make the return matrix
                colnames(est.pars) <- c("est","CIlow", "CIhigh")
                rownames(est.pars) <- c(par1.name,par2.name,paste0("p", 100*ptiles.appended))

                #making the matrix with the actual estimates.
                est.pars[1,] <- quantile(untrans.mcmcs[,1], c(0.5,0.025,0.975))
                est.pars[2,] <- quantile(untrans.mcmcs[,2], c(0.5,0.025,0.975))
                cis.ptiles <- t(apply(mcmc.quantiles,1,function(x) quantile(x,c(0.5,.025,.975))))
                est.pars[3:nrow(est.pars),1:3] <- cis.ptiles

                ## finally get tbhe log-likelihood evaluated at the mean posterior for each parameter                            
                ll <- -loglikhd(pars=dist.optim.transform(dist=dist,est.pars[1:2,1]),dat=data.frame(dat),dist=dist)                                

                rc <- new("cd.fit.mcmc",
                          ests=round(est.pars,3),
                          conv = numeric(),
                          MSG = "",
                          loglik=ll,
                          samples = data.frame(untrans.mcmcs),
                          data=data.frame(dat),
                          dist=dist,
                          inv.hessian = matrix(),
                          est.method = "MCMC",
                          ci.method = "MCMC"
                )

                return(rc)

        } else {
            if (msg != "") msg <- paste0("Error: ",msg)
            stop(sprintf("\n%s\nTry adjusting the starting parameters (init.pars), or changing the optimization method (opt.method).",msg))
        }
    }



##' posterior log likelihood function to pass to MCMCpack sampler
##'
##' @param pars the parameters to calculate the ll at
##' @param dat the date to base it on
##' @param prior.par1 first parameter of each prior
##' @param prior.par2 second parameter of each prior
##' @param dist the distribution the likelihood is being calculated for  
##' 
##' @return the posterior log likelihood
mcmcpack.ll <- function(pars,
                     dat,
                     prior.par1,
                     prior.par2,
                     dist) {
  
  
  
  ## get parameters on untransformed scale
  pars.untrans <- dist.optim.untransform(dist,pars)
  
  
  
  if (dist == "L"){
    ## default gamma on scale param and (inproper) uniform on location
    ll <- tryCatch(-loglikhd(pars,dat,dist) +
                     ## dgamma(pars.untrans[2],shape=par.prior.param1[2],
                     ## rate=par.prior.param2[2],log=T),
                     sum(dnorm(pars,
                               prior.par1,
                               prior.par2,log=T)),
                   error=function(e) {
                     warning("Loglik failure, returning -Inf")
                     return(-Inf)
                   })
    
  } else if (dist == "W"){
    ## using normal prior on the log-shape param and gamma on scale
    ll <- tryCatch(
      -loglikhd(pars,dat,dist) +
        dnorm(pars[1],
              prior.par1[1],
              prior.par2[1],log=T) +
        dgamma(pars.untrans[2],
               shape=prior.par1[2],
               rate=prior.par2[2],log=T),
      error=function(e) {
        warning("Loglik failure, returning -Inf")
        return(-Inf)
      })
  } else if (dist == "G"){
    ## using Jeffery's prior
    ll <- tryCatch(-loglikhd(pars,dat,dist)
                   +
                     log(1/pars.untrans[2]*sqrt(pars.untrans[1]*trigamma(pars.untrans[1])-1))
                   ,
                   error=function(e) {
                     warning("Loglik failure, returning -Inf")
                     return(-Inf)
                   })
    
  } else {
    stop("Sorry, unknown distribution type. Check the 'dist' option.")
  }
  return(ll)
}

##' Does a metropolis hastings for the Erlang distribution 
##' 
##' @param dat the data to fit
##' @param prior.par1 mean of priors. A negative binomial (for shape) and a normal for log(scale)
##' @param prior.par2 dispersion parameters for priors, dispersion for negative binomial, log scale sd for normal
##' @param init.pars the starting parameters on the reporting scale
##' @param verbose how often to print an update
##' @param burnin how many burnin iterations to do
##' @param n.samples the number of samples to keep and report back
##' @param sds the standard deviations for the proposal distribution
##' 
##' @return a matrix of n.samples X 2 parameters, on the estimation scale
##' 
mcmc.erlang <- function (dat, prior.par1, prior.par2, 
                         init.pars, verbose, burnin, n.samples, sds=c(1,1)) {

  #make a logging return matrix
  states <- matrix(nrow=burnin+n.samples, ncol=4)
  colnames(states) <- c("shape", "scale", "ll", "accept")
  
  #calculate the LL for the initial parameters. 
  shape.cur <- init.pars[1]
  scale.cur <- log(init.pars[2])
  ll.cur <- -loglikhd(c(log(shape.cur),scale.cur), dat, dist="G") +
      dnbinom(shape.cur, mu=prior.par1[1], size=prior.par2[1], log=T) +
        dnorm(scale.cur, prior.par1[2], prior.par2[2], log=T)
  
  if(ll.cur==-Inf) {
    stop("zero starting likelihood. check priors are appropriate for Erlang distribution")
  }
  
  states[1,] <- c(shape.cur, scale.cur, ll.cur, 1)
  
  for (i in 2:(burnin+n.samples)) {
    #propose new parameters
    shape.prop <- shape.cur + round(rnorm(1,0,sds[1])) #discretized normal
    scale.prop <- scale.cur + rnorm(1,0,sds[2])
    
    #calculate the liklihood using the new parameters...go ahead and do -Inf if shape is -
    if (shape.prop<0) {
      ll.prop <- -Inf
    } else {
      ll.prop <- -loglikhd(c(log(shape.prop),scale.prop), dat, dist="G") +
        dnbinom(shape.prop, mu=prior.par1[1], size=prior.par2[1], log=T) +
        dnorm(scale.prop, prior.par1[2], prior.par2[2], log=T)
    }
  
    #accept or reject
    l.alpha <- log(runif(1))
    
    
    if ((ll.prop-ll.cur)>l.alpha) {
      #accept
      shape.cur <- shape.prop
      scale.cur <- scale.prop
      ll.cur <- ll.prop
      states[i,] <- c(shape.cur, scale.cur, ll.cur, 1)
    } else {
      #reject
      states[i,] <- c(shape.cur, scale.cur, ll.cur, 0)
    }
  
    if (verbose !=0 & i%%verbose == 0) {
        cat("Erlang MCMC iteration ",i," of ", burnin+n.samples,"\n")
        cat("LL = ", ll.cur, "\n")
        cat("theta = ", c(shape.cur, scale.cur),"\n")
        cat("acceptance rate = ",sum(states[1:i,4])/i,"\n")
    }
  
  }

  rc <- states[(burnin+1):(burnin+n.samples),1:2]
  return(rc)
}
































