dmdirichlet=function(x, mAlpha, mixtureCoef) {
    d=0;
    for (i in 1:nrow(mAlpha)) {
        d=d + mixtureCoef[i]*ddirichlet(x, mAlpha[i,]);
    }
    d
}

# copied from MCMCpack
ddirichlet=function (x, alpha)
{
    dirichlet1 <- function(x, alpha) {
        logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        s <- sum((alpha - 1) * log(x))
        exp(sum(s) - logD)
    }
    if (!is.matrix(x))
        if (is.data.frame(x))
            x <- as.matrix(x)
        else x <- t(x)
    if (!is.matrix(alpha))
        alpha <- matrix(alpha, ncol=length(alpha), nrow=nrow(x),
            byrow=TRUE)
    if (any(dim(x) != dim(alpha)))
        stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
    pd <- vector(length=nrow(x))
    for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,
        ])
    pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
    return(pd)
}


# copied from MCMCpack
rdirichlet=function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}


rmdirichlet=function(mAlpha, mixtureCoef) {
    p=runif(1)
    k=min(which (p<cumsum (mixtureCoef)))
    # above implementation is similar to this, but use rng differently: k=which(rmultinom(1, 1, mixtureCoef)==1)
    rdirichlet(1, mAlpha[k,])
}

my.lbeta=function (x) {
    sum(lgamma(x)) - lgamma(sum(x))
}

modifyDirichlet=function(prior, y) {
    alpha=prior$alpha
    mix.coef=prior$mix.coef
    new.coef=sapply(1:nrow(alpha), function(d){
        alpha0=alpha[d,]
        log(mix.coef[d]) + my.lbeta(alpha0+y) - my.lbeta(alpha0)
    })    
    list (alpha=matrix(t(alpha)+y, byrow=T, ncol=20), 
        mix.coef=exp(new.coef-logSumExp(new.coef))  )
}

# y is aa counts from one column
logIntegrateMixDirichlet=function(y, prior, tau=1) {
    alpha=prior$alpha
    mix.coef=prior$mix.coef
    new.coef=sapply(1:nrow(alpha), function(d){
        alpha=tau*alpha[d,]
        log(mix.coef[d]) + my.lbeta(alpha+y) - my.lbeta(alpha)
    })    
    logSumExp(new.coef)
}

logIntegrateDirichlet=function (y, alpha) {
    my.lbeta(alpha+y) - my.lbeta(alpha)
}
