LambdaOptim <-
function(nsim,bvec,f1,f2,p1,p2,init){

    p <- rdirichlet(nsim,alpha=bvec)	
    pmin <- apply(p,1,min)
    fmin <- -pmin/(1-pmin)	

    lambda0 <- rnorm(nsim)
    
    LambdaPriorChoice <- function(x) {
        lambda <- lambda0 * exp(x[2]) + x[1]
        f <- (exp(lambda)+fmin)/(exp(lambda)+1)
        (quantile(f,probs=p1)-f1)^2+(quantile(f,probs=p2)-f2)^2
    }

    opt <- optim(par=init,fn=LambdaPriorChoice)
    lambdamu <- opt$par[1]
    lambdasd = exp(opt$par[2])
    cat("lambda mu and lambda sd = ",lambdamu,lambdasd,"\n")
    list(lambdamu=lambdamu,lambdasd=lambdasd)
}

