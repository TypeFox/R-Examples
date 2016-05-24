## Make term for sscopu
mkterm.copu <- function(dm,order,symmetry,exclude)
{
    phi0 <- function(x) x-.5
    rk0 <- function(x,y) {
        k2 <- function(x) ((x-.5)^2-1/12)/2
        k4 <- function(x) ((x-.5)^4-(x-.5)^2/2+7/240)/24
        k2(x)*k2(y)-k4(abs(x-y))
    }
    combo <- function(x,order) {
        d <- length(x)
        if (d==order) return(matrix(x,nrow=1))
        if (order==1) z <- matrix(x,ncol=1)
        else {
            wk <- NULL
            for (i in 2:(d-order+2))
                wk <- rbind(wk,combo(x[i:d],order-1))
            z <- cbind(rep(x[1:(d-order+1)],choose((d-1):(order-1),order-1)),wk)
        }
        z
    }
    if (symmetry) {
        permut <- function(x) {
            d <- length(x)
            if (d==1) return(x)
            z <- NULL
            for (i in 1:d) z <- rbind(z,cbind(x[i],permut(x[-i])))
            z
        }
        nphi <- order
        phi <- function(x,nu,env) {
            wk <- env$combo(1:dim(x)[2],nu)
            z <- 0
            for (i in 1:dim(wk)[1]) {
                z.wk <- 1
                for (j in 1:nu) z.wk <- z.wk*env$phi0(x[,wk[i,j]])
                z <- z + z.wk
            }
            z
        }
        perm <- permut(1:dm)
        idx.rk <- cumsum(1:order)
        nrk <- sum(1:order)
        rk <- function(x,y,nu,env,out=TRUE) {
            order <- min((1:length(env$idx.rk))[nu<=env$idx.rk])
            if (order==1) {
                z <- 0
                for (i in 1:dim(env$perm)[1]) {
                    y.wk <- y[,perm[i,]]
                    for (j in 1:dim(x)[2]) {
                        if (out) z <- z + outer(x[,j],y.wk[,j],env$rk0)
                        else z <- z + env$rk0(x[,j],y.wk[,j])
                    }
                }
                return(z)
            }
            wk <- env$combo(1:dim(x)[2],order)
            nu.wk <- nu - env$idx.rk[order-1]
            z <- 0
            for (i in 1:dim(env$perm)[1]) {
                y.wk <- y[,perm[i,]]
                for (j in 1:dim(wk)[1]) {
                    idx.wk <- wk[j,]
                    wk1 <- env$combo(idx.wk,nu.wk)
                    for (k in 1:dim(wk1)[1]) {
                        z.wk <- 1
                        for (ind in idx.wk) {
                            if (ind%in%wk1[k,]) {
                                if (out) z.wk <- z.wk*outer(x[,ind],y.wk[,ind],env$rk0)
                                else z.wk <- z.wk*env$rk0(x[,ind],y.wk[,ind])
                            }
                            else {
                                if (out) z.wk <- z.wk*outer(env$phi0(x[,ind]),
                                                            env$phi0(y.wk[,ind]))
                                else z.wk <- z.wk*env$phi0(x[,ind])*env$phi0(y.wk[,ind])
                            }
                        }
                        z <- z + z.wk
                    }
                }
            }
            z
        }
        env <- list(phi0=phi0,rk0=rk0,combo=combo,perm=perm,idx.rk=idx.rk)
    }
    else {
        idx.phi <- NULL
        nphi <- 0
        for (i in 1:order) {
            wk <- combo(1:dm,i)
            if (!is.null(exclude)) {
                exc.ind <- NULL
                for (j in 1:choose(dm,i)) {
                    exc.wk <- FALSE
                    for (k in 1:dim(exclude)[1])
                        exc.wk <- any(exc.wk,all(exclude[k,]%in%wk[j,]))
                    exc.ind <- c(exc.ind,exc.wk)
                }
                wk <- wk[!exc.ind,,drop=FALSE]
            }
            if (!dim(wk)[1]) next
            for (j in 1:dim(wk)[1]) idx.phi <- rbind(idx.phi,(1:dm)%in%wk[j,])
            nphi <- nphi + dim(wk)[1]
        }
        phi <- function(x,nu,env) {
            z <- 1
            for (i in 1:dm)
                if (env$idx.phi[nu,i]) z <- z*env$phi0(x[,i])
            z
        }
        idx.rk <- NULL
        nrk <- 0
        for (i in 1:order) {
            wk <- combo(1:dm,i)
            if (!is.null(exclude)) {
                exc.ind <- NULL
                for (j in 1:choose(dm,i)) {
                    exc.wk <- FALSE
                    for (k in 1:dim(exclude)[1])
                        exc.wk <- any(exc.wk,all(exclude[k,]%in%wk[j,]))
                    exc.ind <- c(exc.ind,exc.wk)
                }
                wk <- wk[!exc.ind,,drop=FALSE]
            }
            if (!dim(wk)[1]) next
            for (j in 1:dim(wk)[1]) {
                idx.wk <- (1:dm)%in%wk[j,]
                for (k in 1:(2^i-1)) {
                    ind <- k
                    tmp <- NULL
                    for (kk in 1:i) {
                        tmp <- c(tmp,ind%%2)
                        ind <- ind%/%2
                    }
                    wwk <- rep(0,dm)
                    wwk[idx.wk] <- tmp+1
                    idx.rk <- rbind(idx.rk,wwk)
                }
            }
            nrk <- nrk + dim(wk)[1]*(2^i-1)
        }
        rk <- function(x,y,nu,env,out=TRUE) {
            z <- 1
            if (out) {
                for (i in 1:dm) {
                    if (env$idx.rk[nu,i]==1)
                        z <- z*outer(env$phi0(x[,i]),env$phi0(y[,i]))
                    if (env$idx.rk[nu,i]==2)
                        z <- z*outer(x[,i],y[,i],env$rk0)
                }
            }
            else {
                for (i in 1:dm) {
                    if (env$idx.rk[nu,i]==1)
                        z <- z*env$phi0(x[,i])*env$phi0(y[,i])
                    if (env$idx.rk[nu,i]==2) z <- z*env$rk0(x[,i],y[,i])
                }
            }
            z
        }
        env <- list(phi0=phi0,rk0=rk0,idx.phi=idx.phi,idx.rk=idx.rk)
    }
    list(nphi=nphi,phi=phi,nrk=nrk,rk=rk,env=env)
}
