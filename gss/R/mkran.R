## Make random effects for mixed-effect models
mkran <- function(formula,data)
{
    with(data,{
    ## decipher formula
    form.wk <- terms.formula(formula)[[2]]
    if (!("|"%in%strsplit(deparse(form.wk),'')[[1]]))
        stop("gss error in mkran: missing | in grouping formula")
    term.wk <- strsplit(deparse(form.wk),' \\| ')[[1]]
    ## make matrix Z
    z2.wk <- eval(parse(text=term.wk[2]))
    if (!is.factor(z2.wk))
        stop(paste("gss error in mkran: ", term.wk[2], " should be a factor"))
    z <- NULL
    lvl.z2 <- levels(z2.wk)
    for (i in lvl.z2) z <- cbind(z,as.numeric(z2.wk==i))
    ## make sigma function
    if (term.wk[1]=="1") {
        init <- 0
        env <- length(levels(z2.wk))
        fun <- function(zeta,env) diag(10^(-zeta),env)
        sigma <- list(fun=fun,env=env)
    }
    else {
        z1.wk <- eval(parse(text=term.wk[1]))
        if (!is.factor(z1.wk))
            stop(paste("gss error in mkran: ", term.wk[1], " should be a factor"))
        ind <- lvl.wk <- NULL
        nz <- length(lvl.z2)
        nsig <- length(levels(z1.wk))
        for (i in levels(z1.wk)) {
            zz.wk <- z2.wk[z1.wk==i,drop=TRUE]
            ind <- c(ind,list((1:nz)[lvl.z2%in%levels(zz.wk)]))
            lvl.wk <- c(lvl.wk,levels(zz.wk))
        }
        if (max(table(lvl.wk)>1))
            stop("gss error in mkran: ", term.wk[2], " should be nested under ", term.wk[1])
        init <- rep(0, length(levels(z1.wk)))
        env <- list(size=nz,nsig=nsig,ind=ind)
        fun <- function(zeta,env) {
            wk <- rep(0,env$size)
            for (i in 1:env$nsig) wk[env$ind[[i]]] <- 10^(-zeta[i])
            diag(wk)
        }
        sigma <- list(fun=fun,env=env)
    }
    list(z=z,sigma=sigma,init=init)
    })
}
