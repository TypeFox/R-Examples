#Code for the link function in the additive model with exponential fitness
## z=x+y where x, y are fitnesses from Exp(1) distribution
linkfun <- function(z, alpha, beta, gamma=1){
    if(alpha==-1){
        res1<-beta*(gamma/beta)^(1-exp(-z))*(1-log(gamma/beta)*exp(-z));
        return(res1);
    }else{
        tmp1 <-(gamma^(alpha+1) + (beta^(alpha+1) - gamma^(alpha+1))*exp(-z))^(1/(alpha+1));
        tmp2<-(1/(alpha+1))*(gamma^(alpha+1) + (beta^(alpha+1) - gamma^(alpha+1))*exp(-z))^(-alpha/(alpha+1))*(gamma^(alpha+1)-beta^(alpha+1))*exp(-z);
        return(tmp1-tmp2);
    }
}

#' Mean out-degree of a random node  the fitness model
#'
#' Computes the relative mean out-degree of a randomly chosen node
#' given fitness \code{x} in the fitness model implemented in
#' \code{\link{Model.additivelink.exponential.fitness}}. The function
#' returns the mean out-degree divided by n-1.
#'
#' @inheritParams Model.additivelink.exponential.fitness
#' 
#' @export
Model.fitness.meandegree <- function(alpha,beta,gamma=1){
    fsum <-function(x) dgamma(x,shape=2,scale=1)
    integrate(function(x) linkfun(x,alpha=alpha,beta=beta,gamma=gamma)*fsum(x),lower=0,upper=Inf)$value
}


#' Mean out-degree of a node with given fitness in the fitness model
#'
#' Computes the mean out-degree of a node with given fitness \code{x}
#' in the fitness model implemented in
#' \code{\link{Model.additivelink.exponential.fitness}}. The function
#' returns the mean out-degree divided by n-1.
#'
#' @param x Fitness of node. A nonegative number.
#' @inheritParams Model.additivelink.exponential.fitness
#' 
#' @export
Model.fitness.conditionalmeandegree <- function(x,alpha,beta,gamma=1){
    fsum <-function(y) dexp(y)
    integrate(function(y) linkfun(y,alpha=alpha,beta=beta,gamma=gamma)*fsum(y-x),lower=x,upper=Inf)$value
}

#' Prior distribution for eta and zeta in the fitness model
#'
#' Assumes a uniform distribution on the shape parameter \code{zeta} and an
#' exponential distribution on the scale parameter \code{eta}. To be used
#' as prior for \code{\link{Model.additivelink.exponential.fitness}}.
#' 
#' @param shapemin Minimal Value of the shape parameter. Default: 0.75.
#' @param shapemax Maximal Value of the shape parameter. Default: 1.5.
#' @param ratescale  Rate parameter for the prior distribution of the scale parameter. In the model this is on the same scale as the entries of \code{L}
#' @param sdshapeprob Standard deviation for the additivel normally
#' distributed random walk proposal for the shape parameter. Defaults to 0.1.
#' @param sdpropscale Standard deviation for the multiplicative
#' lognormal proposals for the scale parameter.
#'
#' @return list of functions necessary for constructing
#' Metropolis-Hastings updates.
#' @export
Model.fitness.genlambdaparprior <- function(shapemin=0.75,shapemax=1.5,ratescale,sdshapeprob=0.1,sdpropscale=0.1){
    list(r=function()c(shape=runif(1,min=shapemin,max=shapemax),scale=rexp(1,rate=ratescale)),
         logd=function(shape,scale){
             (dunif(shape,min=shapemin,max=shapemax,log=TRUE)+
                  dexp(scale,rate=ratescale,log=TRUE))
         },
         prop=function(shape,scale){
             delta <- rlnorm(1,meanlog=0,sdlog=sdpropscale)
             list(prop=c(shape=rnorm(1,mean=shape,sd=sdshapeprob),
                      scale=scale*delta),
                  logqratio=(dlnorm(1/delta,meanlog=0,sdlog=sdpropscale,log=TRUE)-dlnorm(delta,meanlog=0,sdlog=sdpropscale,log=TRUE)
                             -log(delta)
                             )##correction for Metropolis Hastings
                  )
         })
}

#' Fitness model for liabilities matrix
#'
#' Assumes a diagonal consisting of 0s.
#' 
#' @param n Number of nodes in the model.
#' @param alpha Exponent of the power law of the degree distribution. Must be <0.
#' @param beta Lower endpoint of the relative expected out degree
#' (expected out degree divided by n-1). Must be >=0.
#' @param gamma Upper endpoint of the relative expected out degree
#' (expected out degree divided by n-1). Must be at least beta and at most 1.
#' @param lambdaprior Prior on zeta and eta. For the type of object
#' required see \code{\link{Model.fitness.genlambdaparprior}}.
#' @param sdpropfitness Standard deviation for the log-normally
#' distributed multiplicative proposals for Metropoli-Hastings updates
#' of the fitness. Defaults to \code{1/sqrt{n}}.
#'
#' @return A model to be used by sample_HierarchicalModel. This is a
#' list of functions. It includes a function accrates() that repors acceptance rates for the Metropolis-Hasting steps involved.
#'
#' @examples
#' mod <- Model.additivelink.exponential.fitness(10,alpha=-2.5,beta=0.1,
#'                 lambdaprior=Model.fitness.genlambdaparprior(ratescale=1e3))
#' theta <- mod$rtheta()
#' L <- genL(mod)
#' l <- rowSums(L$L)
#' a <- colSums(L$L)
#' ## increase number of samples and thinning in real examples
#' res <- sample_HierarchicalModel(l=l,a=a,model=mod,nsamples=4,thin=50)
#' mod$accrates()
#' @export
Model.additivelink.exponential.fitness <- function(n,alpha,beta,gamma=1,lambdaprior,sdpropfitness=1/sqrt(n)){
    if (alpha>=0) stop("alpha must be <0.")
    if (beta<=0) stop("beta must be positive.")
    if (beta>gamma) stop("beta must be less or equal to gamma.")
    if (gamma>1) stop("gamma must be <=1")
    if (alpha==-1 && (1-log(gamma/beta)<0))
        stop("Parameter values lead to a negative link function. Need gamma/beta>=exp(-1)")
    if (alpha>-2&&alpha< -1&&gamma/beta<(alpha+2)^(1/(alpha+1)))
        stop("Parameter values lead to a negative link function. Need
gamma/beta>=(alpha+2)^(1/(alpha+1))")
    if (alpha>-1 && gamma/beta>(alpha+2)^(1/(alpha+1)))
        stop("Parameter values lead to a negative link function. Need
gamma/beta<=(alpha+2)^(1/(alpha+1))")
    lf <- linkfun(c(seq(0,10,by=0.01),Inf),alpha=alpha,beta=beta,gamma=gamma)
    if (min(lf)<0)
        stop("Parameter values lead to a negative link function.")
    if (max(lf)<0)
        stop("Parameter values lead to a link probability>1")
    genmatrlambda <- function(theta){
        xtransformed <- qgamma(exp(-theta[ifitness]),
                               shape=theta[ilambdashape],
                               scale=theta[ilambdascale])
        if (any(is.na(xtransformed))) browser()
        outer(xtransformed,xtransformed,FUN="+")
    }
    genmatrp <- function(fitness){
        m <- outer(fitness,fitness,FUN=function(x,y) linkfun(x+y,alpha=alpha,beta=beta,gamma=gamma))
        diag(m) <- 0
        m
    }
    nupdateaccept_fitness <- 0;
    nupdateaccept_lambdapar <- 0
    nupdates <- 0;
    ifitness <- 1:n
    ilambdashape <- n+1
    ilambdascale <- n+2
    ilambdapar <- c(ilambdashape,ilambdascale)
    list(dim=n+2,
         update=function(L,theta){
             loglikelihood <- function(theta){
                 res <- lambdaprior$logd(shape=theta[ilambdashape],scale=theta[ilambdascale])
                 if (res==-Inf) return(res);
                 pakt <- genmatrp(theta[ifitness])
                 lambdaakt <- genmatrlambda(theta)
                 prob <- ifelse(L>0,log(pakt)+log(lambdaakt)-lambdaakt*L,
                                log(1-pakt))
                 diag(prob) <- 0
                 (sum(prob)+sum(dexp(theta[ifitness],log=TRUE))
                  +res
                  )
             }
             nupdates <<- nupdates+1

             thetanew <- theta
             delta <- rlnorm(n,sdlog=sdpropfitness)
             thetanew[ifitness] <- thetanew[ifitness]*delta

             if (all(thetanew[ifitness]>=0)){
                 loglikold <- loglikelihood(theta)
                 if (loglikold==-Inf){
                     if (loglikelihood(thetanew)>-Inf){
                         theta <- thetanew
                         nupdateaccept_fitness <<- nupdateaccept_fitness +1
                     }
                 }else{
                     if (exp(loglikelihood(thetanew)-loglikold
                             +sum(dlnorm(1/delta,sdlog=sdpropfitness,log=TRUE))
                             -sum(dlnorm(delta,sdlog=sdpropfitness,log=TRUE))
                             -sum(log(delta)))
                                 >runif(1)){
                         theta <- thetanew
                         nupdateaccept_fitness <<- nupdateaccept_fitness +1
                     }

                 }
             }
             ## update parameter for lambda
             thetanew <- theta
             lambdaprop <- lambdaprior$prop(shape=theta[ilambdashape],scale=theta[ilambdascale])
             thetanew[ilambdapar] <- lambdaprop$prop
             loglikold <- loglikelihood(theta)
             if (loglikold==-Inf){
                 if (loglikelihood(thetanew)>-Inf){
                     theta <- thetanew
                     nupdateaccept_lambdapar <<- nupdateaccept_lambdapar +1
                 }else{
                     browser()
                 }
             }else{
                 if (exp(loglikelihood(thetanew)-loglikold+lambdaprop$logqratio)
                     >runif(1)){
                     theta <- thetanew
                     nupdateaccept_lambdapar <<- nupdateaccept_lambdapar +1
                 }
             }
             ##
             theta
         },
         matr=function(theta){
             list(p=genmatrp(theta[ifitness]),
                  lambda=genmatrlambda(theta))
         },
         rtheta=function(){
             c(fitness=rexp(n),lambdaprior$r())
         },
         inittheta=function(){
             c(fitness=rep(1,n),lambdaprior$r())
         },
         accrates=function(){
             cat("Acceptance rate for fitness updates:",nupdateaccept_fitness/nupdates,"\n")
             cat("Acceptance rate for lambdapar updates:", nupdateaccept_lambdapar/nupdates,"\n")
         })

}

