## Calculate square error projection from ssden1 objects
project.ssden1 <- function(object,include,drop1=FALSE,...)
{
    ## calculate log(rho) and integrals
    rho1 <- sum(object$rho.int)
    rho2 <- rho1^2-sum(object$rho.int^2)+sum(object$rho.int2)
    ## calculate cross integrals of rho, phi, and rk
    s <- object$int$s
    r <- object$int$r
    s.rho <- object$int$s.rho - s*rho1
    r.rho <- object$int$r.rho - r*rho1
    ## obtain ss, sr, rr
    int2 <- mkint2(object$mf,object$int$var.type,
                   object$id.basis,object$quad,object$terms)
    ss <- int2$ss
    sr <- int2$sr
    rr <- int2$rr
    ## evaluate full model
    d <- object$d
    c <- object$c
    theta <- object$theta
    nq <- length(theta)
    s.eta <- ss%*%d
    r.eta <- tmp <- NULL
    r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
    for (i in 1:nq) {
        tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
        s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
        if (length(d)==1) r.eta.wk <- sr[,,i]*d
        else r.eta.wk <- t(sr[,,i])%*%d
        r.wk <- r.wk + 10^theta[i]*r[,i]
        r.rho.wk <- r.rho.wk + 10^theta[i]*r.rho[,i]
        sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
        for (j in 1:nq) {
            r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
            rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
        }
        r.eta <- cbind(r.eta,r.eta.wk)
    }
    s.eta <- s.eta - s*(sum(s*d)+sum(r.wk*c))
    r.eta <- r.eta - r*(sum(s*d)+sum(r.wk*c))
    ss <- ss - outer(s,s,"*")
    sr.wk <- sr.wk - outer(s,r.wk,"*")
    rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
    rho.eta <- sum(s.rho*d) + sum(r.rho.wk*c)
    eta2 <- sum(c*(rr.wk%*%c)) + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
    mse <- eta2 + rho2-rho1^2 + 2*rho.eta
    ## calculate projection
    rkl <- function(include) {
        inc.wk <- union(names(object$mf),include)
        id.s <- id.q <- NULL
        for (label in inc.wk) {
            if (!any(label==object$terms$labels)) next
            term <- object$terms[[label]]
            if (term$nphi>0) id.s <- c(id.s,term$iphi+(1:term$nphi)-2)
            if (term$nrk>0) id.q <- c(id.q,term$irk+(1:term$nrk)-1)
        }
        ss.wk <- ss[id.s,id.s]
        r.eta.wk <- r.wk <- sr.wk <- rr.wk <- 0
        for (i in id.q) {
            r.eta.wk <- r.eta.wk + 10^theta[i]*r.eta[,i]
            r.wk <- r.wk + 10^theta[i]*r[,i]
            sr.wk <- sr.wk + 10^theta[i]*sr[id.s,,i]
            for (j in id.q) {
                rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
            }
        }
        sr.wk <- sr.wk - outer(s[id.s],r.wk,"*")
        rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
        v <- cbind(rbind(ss.wk,t(sr.wk)),rbind(sr.wk,rr.wk))
        mu <- c(s.eta[id.s],r.eta.wk)
        nn <- length(mu)
        z <- chol(v,pivot=TRUE)
        v <- z
        rkv <- attr(z,"rank")
        m.eps <- .Machine$double.eps
        while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
        if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
        mu <- backsolve(v,mu[attr(z,"pivot")],transpose=TRUE)
        eta2 - sum(mu[1:rkv]^2)
    }
    ## projection
    if (drop1) {
        se <- NULL
        for (i in 1:length(include)) se <- c(se,rkl(include[-i]))
        ratio <- se/mse
        names(se) <- names(ratio) <- include
    }
    else se <- rkl(include)
    ratio <- se/mse
    list(ratio=ratio,se=se)
}

## Calculate integrals of phi and rk for ssden1
mkint2 <- function(mf,type,id.basis,quad,term)
{
    ## Obtain model terms
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt,"variables"))[-1]
    xfacs <- attr(mt,"factors")
    term.labels <- labels(mt)
    vlist <- xvars[as.logical(apply(xfacs,1,sum))]
    ## Create phi and rk
    nbasis <- length(id.basis)
    phi.term <- rk.term <- list(NULL)
    nvar <- length(names(mf))
    ns <- nq <- 0
    for (label in term.labels) {
        ns <- ns+term[[label]]$nphi
        nq <- nq+term[[label]]$nrk
        phi.term[[label]] <- rk.term[[label]] <- list(NULL)
        vlist <- xvars[as.logical(xfacs[,label])]
        x <- mf[,vlist]
        dm <- length(vlist)
        phi <- rk <- NULL
        if (dm==1) {
            type.wk <- type[[vlist]][[1]]
            xx <- mf[id.basis,vlist]
            xmesh <- quad[[vlist]]$pt
            if (type.wk%in%c("nominal","ordinal")) {
                ## factor variable
                if (type.wk=="nominal") fun <- mkrk.nominal(levels(x))
                else fun <- mkrk.ordinal(levels(x))
                if (nlevels(x)>2) {
                    ## rk
                    rk <- fun$fun(xmesh,xx,fun$env,TRUE)
                }
                else {
                    ## phi
                    wk <- as.factor(names(fun$env$code)[1])
                    phi <- fun$fun(xmesh,wk,fun$env)
                }
            }
            if (type.wk=="cubic") {
                ## cubic splines
                range <- type[[vlist]][[2]]
                ## phi
                phi.fun <- mkphi.cubic(range)
                phi <- phi.fun$fun(xmesh,1,phi.fun$env)
                ## rk
                rk.fun <- mkrk.cubic(range)
                rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
            }
            if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
                ## cubic periodic, linear, and linear periodic splines
                range <- type[[vlist]][[2]]
                ## rk
                if (type.wk=="cubic.per") rk.fun <- mkrk.cubic.per(range)
                if (type.wk=="linear") rk.fun <- mkrk.linear(range)
                if (type.wk=="linear.per") rk.fun <- mkrk.linear.per(range)
                if (type.wk=="sphere") rk.fun <- mkrk.sphere(range)
                rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
            }
            if (type.wk=="tp") {
                ## thin-plate splines
                par <- type[[vlist]][[2]]
                order <- par$order
                mesh <- par$mesh
                weight <- par$weight
                if (is.vector(x)) xdim <- 1
                else xdim <- dim(x)[2]
                ## phi
                phi.fun <- mkphi.tp(xdim,order,mesh,weight)
                nphi <- choose(xdim+order-1,xdim)-1
                if (nphi>0) {
                    for (nu in 1:nphi) {
                        phi <- cbind(phi,phi.fun$fun(xmesh,nu,phi.fun$env))
                    }
                }
                ## rk
                rk.fun <- mkrk.tp(xdim,order,mesh,weight)
                rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
            }
            if (type.wk=="custom") {
                ## user-defined
                par <- type[[vlist]][[2]]
                nphi <- par$nphi
                if (nphi>0) {
                    phi.fun <- par$mkphi(par$env)
                    for (nu in 1:nphi) {
                        phi <- cbind(phi,phi.fun$fun(xmesh,nu,phi.fun$env))
                    }
                }
                rk.fun <- par$mkrk(par$env)
                rk <- rk.fun$fun(xmesh,xx,rk.fun$env,TRUE)
            }
            phi.term[[label]][[vlist]] <- phi
            if (is.null(rk)) rk.term[[label]][[vlist]] <- rk
            else {
                nmesh <- length(quad[[vlist]]$wt)
                rk.term[[label]][[vlist]] <- array(rk,c(nmesh,nbasis,1))
            }
        }
        else {
            bin.fac <- n.phi <- phi.list <- rk.list <- NULL
            for (i in 1:dm) {
                type.wk <- type[[vlist[i]]][[1]]
                if (type.wk%in%c("nominal","ordinal")) {
                    ## factor variable
                    if (type.wk=="nominal")
                        rk.wk <- mkrk.nominal(levels(x[[i]]))
                    else rk.wk <- mkrk.ordinal(levels(x[[i]]))
                    phi.wk <- rk.wk
                    n.phi <- c(n.phi,0)
                    bin.fac <- c(bin.fac,!(nlevels(x[[i]])>2))
                }
                if (type.wk=="cubic") {
                    ## cubic or linear splines
                    range <- type[[vlist[i]]][[2]]
                    ## phi
                    phi.wk <- mkphi.cubic(range)
                    n.phi <- c(n.phi,1)
                    ## rk
                    rk.wk <- mkrk.cubic(range)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
                    ## cubic periodic, linear, or linear periodic splines
                    range <- type[[vlist[i]]][[2]]
                    n.phi <- c(n.phi,0)
                    phi.wk <- NULL
                    if (type.wk=="cubic.per") rk.wk <- mkrk.cubic.per(range)
                    if (type.wk=="linear") rk.wk <- mkrk.linear(range)
                    if (type.wk=="linear.per") rk.wk <- mkrk.linear.per(range)
                    if (type.wk=="sphere") rk.wk <- mkrk.sphere(range)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk=="tp") {
                    ## thin-plate splines
                    par <- type[[vlist[i]]][[2]]
                    order <- par$order
                    mesh <- par$mesh
                    weight <- par$weight
                    if (is.vector(x[[i]])) xdim <- 1
                    else xdim <- dim(x[[i]])[2]
                    phi.wk <- mkphi.tp(xdim,order,mesh,weight)
                    n.phi <- c(n.phi,choose(xdim+order-1,xdim)-1)
                    rk.wk <- mkrk.tp(xdim,order,mesh,weight)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk=="custom") {
                    ## user-defined
                    par <- type[[vlist[i]]][[2]]
                    n.phi <- c(n.phi,par$nphi)
                    if (par$nphi>0) phi.wk <- par$mkphi(par$env)
                    else phi.wk <- NULL
                    rk.wk <- par$mkrk(par$env)
                    bin.fac <- c(bin.fac,0)
                }
                phi.list <- c(phi.list,list(phi.wk))
                rk.list <- c(rk.list,list(rk.wk))
            }
            ## phi
            id0 <- names(mf)%in%vlist
            nphi <- term[[label]]$nphi
            iphi <- term[[label]]$iphi
            if (nphi>0) {
                for (nu in 1:nphi) {
                    ind <- nu - 1
                    for (i in 1:dm) {
                        phi.wk <- phi.list[[i]]
                        xmesh <- quad[[vlist[i]]]$pt
                        if (bin.fac[i]) {
                            wk <- as.factor(names(phi.wk$env$code)[1])
                            phi <- phi.wk$fun(xmesh,wk,phi.wk$env)
                        }
                        else {
                            code <- ind%%n.phi[i] + 1
                            ind <- ind%/%n.phi[i]
                            phi <- phi.wk$fun(xmesh,code,phi.wk$env)
                        }
                        phi.term[[label]][[vlist[i]]] <-
                            cbind(phi.term[[label]][[vlist[i]]],phi)
                    }
                }
            }
            ## rk
            n.rk <- ifelse(n.phi,2,1)
            nrk <- prod(n.rk) - as.logical(nphi)
            if (nrk>0) {
                for (nu in 1:nrk) {
                    ind <- nu - !nphi
                    for (i in 1:dm) {
                        code <- ind%%n.rk[i] + 1
                        ind <- ind%/%n.rk[i]
                        xx <- mf[id.basis,vlist[[i]]]
                        xmesh <- quad[[vlist[i]]]$pt
                        if (code==n.rk[i]) {
                            rk.wk <- rk.list[[i]]
                            rk <- rk.wk$fun(xmesh,xx,rk.wk$env,TRUE)
                        }
                        else {
                            rk <- 0
                            phi.wk <- phi.list[[i]]
                            for (j in 1:n.phi[i]) {
                                phix <- phi.wk$fun(xmesh,j,phi.wk$env)
                                phiy <- phi.wk$fun(xx,j,phi.wk$env)
                                rk <- rk + outer(phix,phiy)
                            }
                        }
                        nmesh <- length(quad[[vlist[i]]]$wt)
                        rk.term[[label]][[vlist[i]]] <-
                            array(c(rk.term[[label]][[vlist[i]]],rk),
                                  c(nmesh,nbasis,nu))
                    }
                }
            }
        }
    }
    ## create arrays
    ss <- matrix(1,ns,ns)
    sr <- array(1,c(ns,nbasis,nq))
    rr <- array(1,c(nbasis,nbasis,nq,nq))
    for (label1 in term.labels) {
        if (!term[[label1]]$nphi) id.s1 <- NULL
        else id.s1 <- term[[label1]]$iphi+(1:term[[label1]]$nphi)-2
        if (!term[[label1]]$nrk) id.r1 <- NULL
        else id.r1 <- term[[label1]]$irk+(1:term[[label1]]$nrk)-1
        irk1 <- term[[label1]]$irk
        for (label2 in term.labels) {
            if (!term[[label2]]$nphi) id.s2 <- NULL
            else id.s2 <- term[[label2]]$iphi+(1:term[[label2]]$nphi)-2
            if (!term[[label2]]$nrk) id.r2 <- NULL
            else id.r2 <- term[[label2]]$irk+(1:term[[label2]]$nrk)-1
            irk2 <- term[[label2]]$irk
            for (xlab in names(mf)) {
                wmesh <- quad[[xlab]]$wt
                phi1 <- phi.term[[label1]][[xlab]]
                phi2 <- phi.term[[label2]][[xlab]]
                rk1 <- rk.term[[label1]][[xlab]]
                rk2 <- rk.term[[label2]][[xlab]]
                ## ss
                if (!is.null(id.s1)&!is.null(id.s2)) {
                    if ((!is.null(phi1))&(!is.null(phi2))) {
                        ss[id.s1,id.s2] <- ss[id.s1,id.s2]*(t(wmesh*phi1)%*%phi2)
                    }
                    else {
                        if (!is.null(phi1)) {
                            ss[id.s1,id.s2] <- ss[id.s1,id.s2]*apply(wmesh*matrix(phi1),2,sum)
                        }
                        else {
                            if (!is.null(phi2)) {
                                ss[id.s1,id.s2] <- t(t(ss[id.s1,id.s2])*
                                                     apply(wmesh*matrix(phi2),2,sum))
                            }
                        }
                    }
                }
                ## sr
                if (!is.null(id.s1)&!is.null(id.r2)) {
                    if ((!is.null(phi1))&(!is.null(rk2))) {
                        for (i in id.r2) {
                            sr[id.s1,,i] <- sr[id.s1,,i]*(t(wmesh*phi1)%*%rk2[,,i-irk2+1])
                        }
                    }
                    else {
                        if (!is.null(phi1)) {
                            sr[id.s1,,id.r2] <- sr[id.s1,,id.r2]*apply(wmesh*matrix(phi1),2,sum)
                        }
                        else {
                            if (!is.null(rk2)) {
                                for (i in id.r2) {
                                    sr[id.s1,,i] <- t(t(sr[id.s1,,i])*
                                                      apply(wmesh*rk2[,,i-irk2+1],2,sum))
                                }
                            }
                        }
                    }
                }
                ## rr
                if (!is.null(id.r1)&!is.null(id.r2)) {
                    if ((!is.null(rk1))&(!is.null(rk2))) {
                        for (i in id.r1) {
                            for (j in id.r2) {
                                rr[,,i,j] <- rr[,,i,j]*(t(wmesh*rk1[,,i-irk1+1])%*%rk2[,,j-irk2+1])
                            }
                        }
                    }
                    else {
                        if (!is.null(rk1)) {
                            for (i in id.r1) {
                                rr[,,i,id.r2] <- rr[,,i,id.r2]*apply(wmesh*rk1[,,i-irk1+1],2,sum)
                            }
                        }
                        else {
                            if (!is.null(rk2)) {
                                for (j in id.r2) {
                                    rr[,,id.r1,j] <-
                                        aperm(aperm(rr[,,id.r1,j,drop=FALSE],c(2,1,3,4))*
                                              apply(wmesh*rk2[,,j-irk2+1],2,sum),c(2,1,3,4))
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    list(ss=ss,sr=sr,rr=rr)
}
