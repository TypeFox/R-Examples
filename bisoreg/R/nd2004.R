nd2004 <- function(x, y, M, n.sim=5000, n.thin=1, n.burn=1000){
    n <- length(y)
    ## gamma <- seq(0, 1, l=M + 1)
    gamma <- seq(min(x), max(x), l=M + 1)
    W <- t(sapply(x,
                  function(r)
                      ifelse(r>=gamma[-(M + 1)], pmin(r, gamma[-1]) - gamma[-(M + 1)], 0)
                  ))
    data <- list(
        W=structure(.Data=as.numeric(W),.Dim=c(n,M)),
        y=y, n=n, M=M)
    params <- c("beta0", "beta", "sigsq", "tausq", "delta")
    initfun <- function(){
        list(beta0=rnorm(1, 0, 1),
             betastar=rnorm(M),
             isigsq=rgamma(1, 1, 1),
             delta=rgamma(1, 1, 1))
    }
    modfile <- system.file("bugs", "nd2004.bug", package="bisoreg")
    n.iter <- n.burn + n.thin*n.sim
    bugsout <- openbugs(data=data,
                        inits=initfun,
                        parameters.to.save=params,
                        model.file=modfile,
                        n.chains=3,
                        n.thin=n.thin,
                        n.burnin=n.burn,
                        n.iter=n.iter)
    return(bugsout)
}

