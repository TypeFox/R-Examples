HWEImportSamp <-
function(nsim,nvec,ischoice,lambdamu,lambdasd,alpha,
         gmu=rep(0,length(alpha)),
         gsigma=diag(0,nrow=length(alpha),
         ncol=length(alpha)))
{
    ##priorint <- varprior <- 0
    k <- length(alpha)
    if (length(nvec) != k*(k+1)/2) {
        stop("length mismatch between alpha and nvec")
    }
    liknorm <- lfactorial(sum(nvec)) - sum(lfactorial(nvec))
    if (ischoice==1) {
        PrnH1 <- varterm1 <- 0
        if(gsigma[1,1]==0) {
            stop("HWImportSamp: You need to supply gmu and gsigma")
        }
        ## we simulate for the baseline logits and lambda
        phisamp <- rmvnorm(n=nsim, mean=gmu, sigma=gsigma)
        for (i in 1:nsim){
            pval <- invbaselogit(phisamp[i,-k])$probs
            lambda <- phisamp[i,k]
            pmin <- min(pval)
            fmin <- -pmin/(1-pmin)
            f <- (exp(lambda)+fmin)/(exp(lambda)+1)
            likterm <- MultLogLikP(pval, f, nvec) + liknorm
            ##
            ## Log of the determinant of the (k-1)x(k-1) Jacobean, derivs are:
            ## partial p_1/partial phi_{1},...partial p_1/partial phi_{k-1}
            ## ........
            ## partial p_{k-1}/partial phi1,...partial p_{k-1}/partial phi_{k-1}
            ##
            jac <- diag(pval[-k]) - outer(pval[-k], pval[-k])
            ljack <- log(det(jac))
            ##
            ## NB We do not need to calculate a Jacobian term for
            ## lambda = phisamp[i,k] as this is generated on the correct
            ## scale for the prior.
            ##
            prterm1 <- log(ddirichlet(pval, alpha=alpha)) + ljack
            prterm2 <- dnorm(lambda, mean=lambdamu, sd=lambdasd, log=TRUE)
            gterm <- dmvnorm(phisamp[i,1:k],mean=gmu,sigma=gsigma,log=TRUE)
            expterm <- exp(likterm+prterm1+prterm2-gterm)
            ## expprior <- exp(prterm1+prterm2-gterm)
            PrnH1 <- PrnH1 + expterm
            ## priorint <- priorint + expprior
            varterm1 <- varterm1 + expterm^2
            ## varprior <- varprior + expprior^2
            if (i %% 1000 == 0) cat("Samples = ",i,"\n")
        }
        ## priorint <- priorint/nsim
        ## varprior <- (varprior/nsim - priorint^2)/nsim
        ##cat("nsim prior constant (se) 95% interval = ",nsim,priorint,"(",sqrt(varprior),")",priorint-1.96*sqrt(varprior),priorint+1.96*sqrt(varprior),"\n")
    }
    if (ischoice==2){

        pval <- rdirichlet(nsim, alpha=alpha)
        lambdaval <- rnorm(nsim, mean=lambdamu, sd=lambdasd)

        minp <- apply(pval, 1, min)
        minf <- -minp/(1-minp)
        f <- (exp(lambdaval)+minf)/(exp(lambdaval)+1)
        X <- cbind(f, pval)
        LLFUN <- function(x) MultLogLikP(x[-1], x[1], nvec)
        likterm <- apply(X, 1, LLFUN) + liknorm
        expterm <- exp(likterm)

        PrnH1 <- sum(expterm)
        varterm1 <- sum(expterm^2)
    }
    PrnH1 <- PrnH1/nsim
    varest <- (varterm1/nsim - PrnH1^2)/nsim
    cat("nsim norm constant (se) 95% interval:\n")
    cat(nsim,PrnH1,"(",sqrt(varest),")",PrnH1-1.96*sqrt(varest),PrnH1+1.96*sqrt(varest),"\n")
    list(PrnH1=PrnH1,varest=varest)
}
