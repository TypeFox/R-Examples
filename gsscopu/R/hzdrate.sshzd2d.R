hzdrate.sshzd2d <- ## Compute hazard rate
function(object,time,covariates=NULL) {
    if (is.vector(time)) time <- rbind(NULL,time)
    s1 <- survexp.sshzd2d(object,time[,1],covariates,1)
    s2 <- survexp.sshzd2d(object,time[,2],covariates,2)
    s12 <- survexp.sshzd2d(object,time,covariates)
    wk <- data.frame(time[,1],time[,2])
    names(wk) <- c(object$hzd1$tname,object$hzd2$tname)
    wk <- cbind(wk,covariates)
    h1 <- hzdrate.sshzd(object$hzd1,wk)
    h2 <- hzdrate.sshzd(object$hzd2,wk)
    as.vector(s1*h1*s2*h2*dsscopu(object$copu,cbind(s1,s2))/s12)
}

survexp.sshzd2d <- ## Compute survival function
function(object,time,covariates=NULL,job=3) {
    if (!(job%in%1:3)) stop("gss error in survexp.sshzd2d: job must be 1, 2, or 3")
    if (is.vector(time)&(job==3)) time <- rbind(NULL,time)
    if (job!=3) {
        if (!is.null(covariates)) {
            nt <- dim(covariates)[1]
            if (nt==1)
                z <- switch(job,survexp.sshzd(object$hzd1,time,covariates),
                            survexp.sshzd(object$hzd2,time,covariates))
            else {
                if (nt!=length(time))
                    stop("gss error in survexp.sshzd2d: time and covariates must match in size")
                z <- NULL
                for (i in 1:nt)
                    z <- c(z,switch(job,
                                    survexp.sshzd(object$hzd1,time[i],
                                                  covariates[i,,drop=FALSE]),
                                    survexp.sshzd(object$hzd2,time[i],
                                                  covariates[i,,drop=FALSE])))
            }
        }
        else z <- switch(job,survexp.sshzd(object$hzd1,time),
                         survexp.sshzd(object$hzd2,time))
    }
    else {
        ## Set up quadrature
        hsz <- 40
        qdsz <- 2*hsz
        qd <- gauss.quad(qdsz,c(0,1))
        gap <- diff(qd$pt)
        g.wk <- gap[hsz]/2
        for (i in 1:(hsz-2)) g.wk <- c(g.wk,gap[hsz+i]-g.wk[i])
        g.wk <- 2*g.wk
        g.wk <- c(g.wk,1/2-sum(g.wk))
        gap[hsz:1] <- gap[hsz+(1:hsz)] <- g.wk
        brk <- cumsum(c(0,gap))[-(qdsz+1)]
        qd.pt <- cbind(rep(qd$pt,qdsz),rep(qd$pt,rep(qdsz,qdsz)))
        d.qd <- matrix(dsscopu(object$copu,qd.pt),qdsz,qdsz)
        if (!is.null(covariates)) {
            nt <- dim(covariates)[1]
            if (nt==1) {
                s1 <- survexp.sshzd(object$hzd1,time[,1],covariates)
                s2 <- survexp.sshzd(object$hzd2,time[,2],covariates)
                z <- NULL
                for (i in 1:dim(time)[1]) {
                    ind1 <- (1:qdsz)[brk<s1[i]]
                    id.mx1 <- max(ind1)
                    wt1 <- qd$wt[ind1]
                    wt1[id.mx1] <- wt1[id.mx1]*(s1[i]-brk[id.mx1])/gap[id.mx1]
                    ind2 <- (1:qdsz)[brk<s2[i]]
                    id.mx2 <- max(ind2)
                    wt2 <- qd$wt[ind2]
                    wt2[id.mx2] <- wt2[id.mx2]*(s2[i]-brk[id.mx2])/gap[id.mx2]
                    z <- c(z,sum(d.qd[ind1,ind2]*outer(wt1,wt2)))
                }
            }
            else {
                if (nt!=dim(time)[1])
                    stop("gss error in survexp.sshzd2d: time and covariates must match in size")
                z <- NULL
                for (i in 1:nt) {
                    s1 <- survexp.sshzd(object$hzd1,time[i,1],covariates[i,,drop=FALSE])
                    ind1 <- (1:qdsz)[brk<s1]
                    id.mx1 <- max(ind1)
                    wt1 <- qd$wt[ind1]
                    wt1[id.mx1] <- wt1[id.mx1]*(s1-brk[id.mx1])/gap[id.mx1]
                    s2 <- survexp.sshzd(object$hzd2,time[i,1],covariates[i,,drop=FALSE])
                    ind2 <- (1:qdsz)[brk<s2]
                    id.mx2 <- max(ind2)
                    wt2 <- qd$wt[ind2]
                    wt2[id.mx2] <- wt2[id.mx2]*(s2-brk[id.mx2])/gap[id.mx2]
                    z <- c(z,sum(d.qd[ind1,ind2]*outer(wt1,wt2)))
                }
            }
        }
        else {
            s1 <- survexp.sshzd(object$hzd1,time[,1])
            s2 <- survexp.sshzd(object$hzd2,time[,2])
            z <- NULL
            for (i in 1:dim(time)[1]) {
                ind1 <- (1:qdsz)[brk<s1[i]]
                id.mx1 <- max(ind1)
                wt1 <- qd$wt[ind1]
                wt1[id.mx1] <- wt1[id.mx1]*(s1[i]-brk[id.mx1])/gap[id.mx1]
                ind2 <- (1:qdsz)[brk<s2[i]]
                id.mx2 <- max(ind2)
                wt2 <- qd$wt[ind2]
                wt2[id.mx2] <- wt2[id.mx2]*(s2[i]-brk[id.mx2])/gap[id.mx2]
                z <- c(z,sum(d.qd[ind1,ind2]*outer(wt1,wt2)))
            }
        }
    }
    z
}
