## Make RK for thin-plate splines
mkrk.tp <- function(dm,order,mesh,weight=1)
{
    ## Check inputs
    if (!((2*order>dm)&(dm>=1))) {
        stop("gss error: thin-plate spline undefined for the parameters")
    }
    if (xor(is.vector(mesh),dm==1)
        |xor(is.matrix(mesh),dm>=2)) {
        stop("gss error in mkrk.tp: mismatched inputs")
    }
    if ((min(weight)<0)|(max(weight)<=0)) {
        stop("gss error in mkrk.tp: negative weights")
    }
    ## Set weights
    if (is.vector(mesh)) N <- length(mesh)
    else N <- dim(mesh)[1]
    weight <- rep(weight,len=N)
    weight <- sqrt(weight/sum(weight))
    ## Obtain orthonormal basis
    phi.p <- mkphi.tp.p(dm,order)
    nnull <- choose(dm+order-1,dm)
    s <- NULL
    for (nu in 1:nnull) s <- cbind(s,phi.p$fun(mesh,nu,phi.p$env))
    s <- qr(weight*s)
    if (s$rank<nnull) {
        stop("gss error in mkrk.tp: insufficient normalizing mesh for thin-plate spline")
    }
    q <- qr.Q(s)
    r <- qr.R(s)
    ## Set Q^{T}E(|u_{i}-u_{j}|)Q
    rk.p <- mkrk.tp.p(dm,order)
    pep <- weight*t(weight*rk.p$fun(mesh,mesh,rk.p$env,out=TRUE))
    pep <- t(q)%*%pep%*%q
    ## Create the environment
    env <- list(dim=dm,order=order,weight=weight,
                phi.p=phi.p,rk.p=rk.p,q=q,r=r,mesh=mesh,pep=pep)
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ## Check inputs
        if (env$dim==1) {
            if (!(is.vector(x)&is.vector(y))) {
                stop("gss error in rk: inputs are of wrong types")
            }
            nx <- length(x)
            ny <- length(y)
        }
        else {
            if (is.vector(x)) x <- t(as.matrix(x))
            if (env$dim!=dim(x)[2]) {
                stop("gss error in rk: inputs are of wrong dimensions")
            }
            nx <- dim(x)[1]
            if (is.vector(y)) y <- t(as.matrix(y))
            if (env$dim!=dim(y)[2]) {
                stop("gss error in rk: inputs are of wrong dimensions")
            }
            ny <- dim(y)[1]
        }
        ## Return the results
        nnull <- choose(env$dim+env$order-1,env$dim)
        if (outer.prod) {
            phix <- phiy <- NULL
            for (nu in 1:nnull) {
                phix <- rbind(phix,env$phi.p$fun(x,nu,env$phi.p$env))
                phiy <- rbind(phiy,env$phi.p$fun(y,nu,env$phi.p$env))
            }
            phix <- backsolve(env$r,phix,transpose=TRUE)
            phiy <- backsolve(env$r,phiy,transpose=TRUE)
            ex <- env$rk.p$fun(env$mesh,x,env$rk.p$env,out=TRUE)
            ex <- env$weight*ex
            ex <- t(env$q)%*%ex
            ey <- env$rk.p$fun(env$mesh,y,env$rk.p$env,out=TRUE)
            ey <- env$weight*ey
            ey <- t(env$q)%*%ey
            env$rk.p$fun(x,y,env$rk.p$env,out=TRUE)-t(phix)%*%ey-
                t(ex)%*%phiy+t(phix)%*%env$pep%*%phiy
        }
        else {
            N <- max(nx,ny)
            phix <- phiy <- NULL
            for (nu in 1:nnull) {
                phix <- rbind(phix,env$phi.p$fun(x,nu,env$phi.p$env))
                phiy <- rbind(phiy,env$phi.p$fun(y,nu,env$phi.p$env))
            }
            phix <- backsolve(env$r,phix,transpose=TRUE)
            phix <- matrix(phix,nnull,N)
            phiy <- backsolve(env$r,phiy,transpose=TRUE)
            phiy <- matrix(phiy,nnull,N)
            ex <- env$rk.p$fun(env$mesh,x,env$rk.p$env,out=TRUE)
            ex <- env$weight*ex
            ex <- t(env$q)%*%ex
            ex <- matrix(ex,nnull,N)
            ey <- env$rk.p$fun(env$mesh,y,env$rk.p$env,out=TRUE)
            ey <- env$weight*ey
            ey <- t(env$q)%*%ey
            ey <- matrix(ey,nnull,N)
            fn1 <- function(x,n) x[1:n]%*%x[n+(1:n)]
            fn2 <- function(x,pep,n) t(x[1:n])%*%pep%*%x[n+(1:n)]
            env$rk.p$fun(x,y,env$rk.p$env)-apply(rbind(phix,ey),2,fn1,nnull)-
                apply(rbind(phiy,ex),2,fn1,nnull)+
                    apply(rbind(phix,phiy),2,fn2,env$pep,nnull)
        }
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make phi function for thin-plate splines
mkphi.tp <-  function(dm,order,mesh,weight)
{
    ## Check inputs
    if (!((2*order>dm)&(dm>=1))) {
        stop("gss error: thin-plate spline undefined for the parameters")
    }
    if (xor(is.vector(mesh),dm==1)
        |xor(is.matrix(mesh),dm>=2)) {
        stop("gss error in mkphi.tp: mismatched inputs")
    }
    if ((min(weight)<0)|(max(weight)<=0)) {
        stop("gss error in mkphi.tp: negative weights")
    }
    ## Set weights
    if (is.vector(mesh)) N <- length(mesh)
    else N <- dim(mesh)[1]
    weight <- rep(weight,len=N)
    weight <- sqrt(weight/sum(weight))
    ## Create the environment
    phi.p <- mkphi.tp.p(dm,order)
    nnull <- choose(dm+order-1,dm)
    s <- NULL
    for (nu in 1:nnull) s <- cbind(s,phi.p$fun(mesh,nu,phi.p$env))
    s <- qr(weight*s)
    if (s$rank<nnull) {
        stop("gss error in mkphi: insufficient normalizing mesh for thin-plate spline")
    }
    r <- qr.R(s)
    env <- list(dim=dm,order=order,phi.p=phi.p,r=r)
    ## Create the phi function
    fun <- function(x,nu,env) {
        nnull <- choose(env$dim+env$order-1,env$dim)
        phix <- NULL
        for(i in 1:nnull)
            phix <- rbind(phix,env$phi.p$fun(x,i,env$phi.p$env))
        t(backsolve(env$r,phix,transpose=TRUE))[,nu+1]
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make pseudo RK for thin-plate splines
mkrk.tp.p <- function(dm,order)
{
    ## Check inputs
    if (!((2*order>dm)&(dm>=1))) {
        stop("gss error: thin-plate spline undefined for the parameters")
    }
    ## Create the environment
    if (dm%%2) {                    
        theta <- gamma(dm/2-order)/2^(2*order)/pi^(dm/2)/gamma(order)
    }
    else {
        theta <- (-1)^(dm/2+order+1)/2^(2*order-1)/pi^(dm/2)/
            gamma(order)/gamma(order-dm/2+1)
    }
    env <- list(dim=dm,order=order,theta=theta)
    ## Create the rk.p function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ## Check inputs
        if (env$dim==1) {
            if (!(is.vector(x)&is.vector(y))) {
                stop("gss error in rk: inputs are of wrong types")
            }
        }
        else {
            if (is.vector(x)) x <- t(as.matrix(x))
            if (env$dim!=dim(x)[2]) {
                stop("gss error in rk: inputs are of wrong dimensions")
            }
            if (is.vector(y)) y <- t(as.matrix(y))
            if (env$dim!=dim(y)[2]) {
                stop("gss error in rk: inputs are of wrong dimensions")
            }
        }
        ## Return the results
        if (outer.prod) {               
            if (env$dim==1) {
                fn1 <- function(x,y) abs(x-y)
                d <- outer(x,y,fn1)
            }
            else {
                fn2 <- function(x,y) sqrt(sum((x-y)^2))
                d <- NULL
                for (i in 1:dim(y)[1]) d <- cbind(d,apply(x,1,fn2,y[i,]))
            }
        }
        else {
            if (env$dim==1) d <- abs(x-y)
            else {
                N <- max(dim(x)[1],dim(y)[1])
                x <- t(matrix(t(x),env$dim,N))
                y <- t(matrix(t(y),env$dim,N))
                fn <- function(x) sqrt(sum(x^2))
                d <- apply(x-y,1,fn)
            }
        }
        power <- 2*env$order-env$dim
        switch(1+env$dim%%2,
               env$theta*d^power*log(ifelse(d>0,d,1)),
               env$theta*d^power)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make pseudo phi function for thin-plate splines
mkphi.tp.p <- function(dm,order)
{
    ## Check inputs
    if (!((2*order>dm)&(dm>=1))) {
        stop("gss error: thin-plate spline undefined for the parameters")
    }
    ## Create the environment
    pol.code <- NULL
    for (i in 0:(order^dm-1)) {
        ind <- i; code <- NULL
        for (j in 1:dm) {
            code <- c(code,ind%%order)
            ind <- ind%/%order
        }
        if (sum(code)<order) pol.code <- cbind(pol.code,code)
    }
    env <- list(dim=dm,pol.code=pol.code)
    ## Create the phi function  
    fun <- function(x,nu,env) {
        if (env$dim==1) x <- as.matrix(x)
        if (env$dim!=dim(x)[2]) {
            stop("gss error in phi: inputs are of wrong dimensions")
        }
        apply(t(x)^env$pol.code[,nu],2,prod)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for spherical splines
mkrk.sphere <- function(order)
{
    ## Create the environment
    env <- list(order=order)
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (is.vector(x)) x <- t(as.matrix(x))
        if (is.vector(y)) y <- t(as.matrix(y))
        if (!(is.matrix(x)&is.matrix(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((dim(x)[2]!=2)|(dim(y)[2]!=2)) {
            stop("gss error in rk: inputs are of wrong dimensions")
        }
        if ((max(abs(x[,1]),abs(y[,1]))>90)|(max(abs(x[,2]),abs(y[,2]))>180)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Convert to radian
        lat.x <- x[,1]/180*pi; lon.x <- x[,2]/180*pi
        lat.y <- y[,1]/180*pi; lon.y <- y[,2]/180*pi
        ##% Return the result
        rk <- function(lat.x,lon.x,lat.y,lon.y,order) {
            z <- cos(lat.x)*cos(lat.y)*cos(lon.x-lon.y)+sin(lat.x)*sin(lat.y)
            W <- ifelse(z<1-10^(-10),(1-z)/2,0)
            A <- ifelse(W>0,log(1+1/sqrt(W)),0)
            C <- ifelse(W>0,2*sqrt(W),0)
            switch(order-1,
                   (A*4*W*(3*W-1)+6*W*(1-C)+1)/2,
                   (W*W*(A*((840*W-720)*W+72)+420*W*(1-C)+220*C-150)-4*W+3)/12,
                   (W*W*W*(A*(((27720*W-37800)*W+12600)*W-600)+
                           (13860*(1-C)*W+14280*C-11970)*W-2772*C+1470)+
                    15*W*W-3*W+5)/30) - 1/(2*order-1)
        }
        if (outer.prod) {
            zz <- NULL
            for (i in 1:length(lat.y))
                zz <- cbind(zz,rk(lat.x,lon.x,lat.y[i],lon.y[i],env$order))
        }
        else zz <- rk(lat.x,lon.x,lat.y,lon.y,env$order)
        zz
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}
