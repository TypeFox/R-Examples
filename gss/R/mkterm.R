## Make phi and rk for model terms
mkterm <- function(mf,type)
{
    ## Obtain model terms
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt,"variables"))[-1]
    xfacs <- attr(mt,"factors")
    term.labels <- labels(mt)
    if (attr(attr(mf,"terms"),"intercept"))
        term.labels <- c("1",term.labels)
    ## Backward compatibility
    vlist <- xvars[as.logical(apply(xfacs,1,sum))]
    if (!is.null(type)&&!is.list(type)&&(type%in%c("cubic","linear","tp"))) {
        type.wk <- type
        type <- NULL
        for (xlab in vlist) type[[xlab]] <- type.wk
    }
    ## Set types for marginals    
    var.type <- NULL
    for (xlab in vlist) {
        x <- mf[,xlab]
        if (!is.null(type[[xlab]])) {
            ## Check consistency and set default parameters
            type.wk <- type[[xlab]][[1]]
            if
            (!(type.wk%in%c("ordinal","nominal","cubic","linear","per","trig",
                            "cubic.per","linear.per","tp","sphere","custom")))
                stop("gss error in mkterm: unknown type")
            if (type.wk%in%c("ordinal","nominal")) {
                par.wk <- NULL
                if (!is.factor(x))
                    stop("gss error in mkterm: wrong type")
            }
            if (type.wk%in%c("cubic","linear")) {
                if (length(type[[xlab]])==1) {
                    mn <- min(x)
                    mx <- max(x)
                    par.wk <- c(mn,mx)+c(-1,1)*.05*(mx-mn)
                }
                else par.wk <- type[[xlab]][[2]]
                if (is.factor(x)|!is.vector(x))
                    stop("gss error in mkterm: wrong type")
            }
            if (type.wk%in%c("per","cubic.per","linear.per","trig")) {
                if (type.wk=="per") type.wk <- "cubic.per"
                if (length(type[[xlab]])==1)
                    stop("gss error in mkterm: missing domain of periodicity")
                else par.wk <- type[[xlab]][[2]]
                if (is.factor(x)|!is.vector(x))
                    stop("gss error in mkterm: wrong type")
            }
            if (type.wk=="tp") {
                if (length(type[[xlab]])==1)
                    par.wk <- list(order=2,mesh=x,weight=1)
                else {
                    par.wk <- par.wk1 <- type[[xlab]][[2]]
                    if (length(par.wk1)==1)
                        par.wk <- list(order=par.wk1,mesh=x,weight=1)
                    if (is.null(par.wk$mesh)) par.wk$mesh <- x
                    if (is.null(par.wk$weight)) par.wk$weight <- 1
                }
                if (dim(as.matrix(x))[2]!=dim(as.matrix(par.wk$mesh))[2])
                    stop("gss error in mkterm: wrong dimension in normalizing mesh")
            }
            if (type.wk=="sphere") {
                if (length(type[[xlab]])==1)
                    par.wk <- 2
                else par.wk <- type[[xlab]][[2]]
                if (!(par.wk%in%(2:4)))
                    stop("gss error in mkterm: spherical order not implemented")
            }
            if (type.wk=="custom") par.wk <- type[[xlab]][[2]]
        }
        else {
            ## Set default types
            if (is.factor(x)) {
                ## categorical variable
                if (is.ordered(x)) type.wk <- "ordinal"
                else type.wk <- "nominal"
                par.wk <- NULL
            }
            else {
                ## numerical variable
                if (is.vector(x)) {
                    type.wk <- "cubic"
                    mn <- min(x)
                    mx <- max(x)
                    par.wk <- c(mn,mx)+c(-1,1)*.05*(mx-mn)
                }
                else {
                    type.wk <- "tp"
                    par.wk <- list(order=2,mesh=x,weight=1)
                }
            }
        }
        var.type[[xlab]] <- list(type.wk,par.wk)
    }
    ## Create the phi and rk functions
    term <- list(labels=term.labels)
    iphi.wk <- 1
    irk.wk <- 1
    for (label in term.labels) {
        iphi <- irk <- phi <- rk <- NULL
        if (label=="1") {
            ## the constant term
            iphi <- iphi.wk
            iphi.wk <- iphi.wk + 1
            term[[label]] <- list(iphi=iphi,nphi=1,nrk=0)
            next
        }
        vlist <- xvars[as.logical(xfacs[,label])]
        x <- mf[,vlist]
        dm <- length(vlist)
        if (dm==1) {
            type.wk <- var.type[[vlist]][[1]]
            if (type.wk%in%c("nominal","ordinal")) {
                ## factor variable
                if (type.wk=="nominal") fun.env <- mkrk.nominal(levels(x))
                else fun.env <- mkrk.ordinal(levels(x))
                if (nlevels(x)>2) {
                    ## phi
                    nphi <- 0
                    ## rk
                    rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                        env$fun(x,y,env$env,outer.prod)
                    }
                    nrk <- 1
                    irk <- irk.wk
                    irk.wk <- irk.wk + nrk
                    rk <- list(fun=rk.fun,env=fun.env)
                }
                else {
                    ## phi
                    phi.fun <- function(x,nu=1,env) {
                        wk <- as.factor(names(env$env$code)[1])
                        env$fun(x,wk,env$env)
                    }
                    nphi <- 1
                    iphi <- iphi.wk
                    iphi.wk <- iphi.wk + nphi
                    phi <- list(fun=phi.fun,env=fun.env)
                    ## rk
                    nrk <- 0
                }
            }
            if (type.wk=="cubic") {
                ## cubic splines
                range <- var.type[[vlist]][[2]]
                ## phi
                phi.env <- mkphi.cubic(range)
                phi.fun <- function(x,nu=1,env) env$fun(x,nu,env$env)
                nphi <- 1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
                ## rk
                rk.env <- mkrk.cubic(range)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
            if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
                ## cubic periodic, linear, and linear periodic splines
                range <- var.type[[vlist]][[2]]
                ## phi
                nphi <- 0
                ## rk
                if (type.wk=="cubic.per") rk.env <- mkrk.cubic.per(range)
                if (type.wk=="linear") rk.env <- mkrk.linear(range)
                if (type.wk=="linear.per") rk.env <- mkrk.linear.per(range)
                if (type.wk=="sphere") rk.env <- mkrk.sphere(range)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
            if (type.wk=="trig") {
                ## trigonometric splines
                range <- var.type[[vlist]][[2]]
                ## phi
                phi.env <- mkphi.trig(range)
                phi.fun <- function(x,nu=1,env) env$fun(x,nu,env$env)
                nphi <- 2
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
                ## rk
                rk.env <- mkrk.trig(range)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
            if (type.wk=="tp") {
                ## thin-plate splines
                par <- var.type[[vlist]][[2]]
                order <- par$order
                mesh <- par$mesh
                weight <- par$weight
                if (is.vector(x)) xdim <- 1
                else xdim <- dim(x)[2]
                ## phi
                phi.env <- mkphi.tp(xdim,order,mesh,weight)
                phi.fun <- function(x,nu=1,env) {
                    env$fun(x,nu,env$env)
                }
                nphi <- choose(xdim+order-1,xdim)-1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
                ## rk
                rk.env <- mkrk.tp(xdim,order,mesh,weight)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
            if (type.wk=="custom") {
                ## user-defined
                par <- var.type[[vlist]][[2]]
                nphi <- par$nphi
                if (nphi>0) {
                    phi.env <- par$mkphi(par$env)
                    phi.fun <- function(x,nu=1,env) {
                        env$fun(x,nu,env$env)
                    }
                    iphi <- iphi.wk
                    iphi.wk <- iphi.wk + nphi
                    phi <- list(fun=phi.fun,env=phi.env)
                }
                rk.env <- par$mkrk(par$env)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
        }
        else {
            bin.fac <- n.phi <- phi.list <- rk.list <- NULL
            for (i in 1:dm) {
                type.wk <- var.type[[vlist[i]]][[1]]
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
                    ## cubic splines
                    range <- var.type[[vlist[i]]][[2]]
                    ## phi
                    phi.wk <- mkphi.cubic(range)
                    n.phi <- c(n.phi,1)
                    ## rk
                    rk.wk <- mkrk.cubic(range)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk%in%c("cubic.per","linear","linear.per","sphere")) {
                    ## cubic periodic, linear, or linear periodic splines
                    range <- var.type[[vlist[i]]][[2]]
                    n.phi <- c(n.phi,0)
                    phi.wk <- NULL
                    if (type.wk=="cubic.per") rk.wk <- mkrk.cubic.per(range)
                    if (type.wk=="linear") rk.wk <- mkrk.linear(range)
                    if (type.wk=="linear.per") rk.wk <- mkrk.linear.per(range)
                    if (type.wk=="sphere") rk.wk <- mkrk.sphere(range)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk=="trig") {
                    ## trigonometric splines
                    range <- var.type[[vlist[i]]][[2]]
                    ## phi
                    phi.wk <- mkphi.trig(range)
                    n.phi <- c(n.phi,2)
                    ## rk
                    rk.wk <- mkrk.trig(range)
                    bin.fac <- c(bin.fac,0)
                }
                if (type.wk=="tp") {
                    ## thin-plate splines
                    par <- var.type[[vlist[i]]][[2]]
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
                    par <- var.type[[vlist[i]]][[2]]
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
            if (!all(as.logical(n.phi+bin.fac))) nphi <- 0
            else {
                phi.env <- list(dim=dm,phi=phi.list,n.phi=n.phi,bin.fac=bin.fac)
                phi.fun <- function(x,nu=1,env) {
                    ind <- nu - 1
                    z <- 1
                    for (i in 1:env$dim) {
                        if (env$bin.fac[i]) {
                            wk <- as.factor(names(env$phi[[i]]$env$code)[1])
                            z <- z * env$phi[[i]]$fun(x[[i]],wk,env$phi[[i]]$env)
                        }
                        else {
                            code <- ind%%env$n.phi[i] + 1
                            ind <- ind%/%env$n.phi[i]
                            z <- z * env$phi[[i]]$fun(x[[i]],code,env$phi[[i]]$env)
                        }
                    }
                    z
                }
                nphi <- prod(n.phi+bin.fac)
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
            }
            ## rk
            rk.env <- list(dim=dm,n.phi=n.phi,nphi=nphi,
                           phi=phi.list,rk=rk.list)
            rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                n.rk <- ifelse(env$n.phi,2,1)
                ind <- nu - !env$nphi
                z <- 1
                for (i in 1:env$dim) {
                    code <- ind%%n.rk[i] + 1
                    ind <- ind%/%n.rk[i]
                    if (code==n.rk[i]) {
                        z <- z * env$rk[[i]]$fun(x[[i]],y[[i]],
                                                 env$rk[[i]]$env,outer.prod)
                    }
                    else {
                        z.wk <- 0
                        for (j in 1:env$n.phi[i]) {
                            phix <- env$phi[[i]]$fun(x[[i]],j,env$phi[[i]]$env)
                            phiy <- env$phi[[i]]$fun(y[[i]],j,env$phi[[i]]$env)
                            if (outer.prod) z.wk <- z.wk + outer(phix,phiy)
                            else z.wk <- z.wk + phix * phiy
                        }
                        z <- z * z.wk
                    }
                }
                z
            }
            n.rk <- ifelse(n.phi,2,1)
            nrk <- prod(n.rk) - as.logical(nphi)
            irk <- irk.wk
            irk.wk <- irk.wk + nrk
            rk <- list(fun=rk.fun,env=rk.env)
        }
        term[[label]] <- list(vlist=vlist,
                              iphi=iphi,nphi=nphi,phi=phi,
                              irk=irk,nrk=nrk,rk=rk)
    }
    term
}
