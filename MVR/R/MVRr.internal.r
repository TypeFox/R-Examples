##########################################################################################################################################
# MVR
##########################################################################################################################################

##########################################################################################################################################
# INTERNAL SUBROUTINES
# (never to be called by end-user)
##########################################################################################################################################

##########################################################################################################################################
################
# Usage         :   is.empty(x)
################
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
################
#
##########################################################################################################################################

is.empty <- function(x) {
    if (is.vector(x)) {
        if((length(x) == 0) || (x == "")) {
            return(T)
        } else {
            return(F)
        }
    } else if (is.matrix(x) || is.data.frame(x)) {
        return( ((nrow(x) == 0) || (ncol(x) == 0)) )
    } else {
        return( ((length(x) == 0) || (x == "")) )
    }
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :   is.valid(x, def, ng)
################
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
################
#
##########################################################################################################################################

is.valid <- function(x, def, ng) {
    ok <- TRUE
    for (g in 1:ng) {
        ok <- (ok && (length(unique(x[def[[g]]])) > 1))
    }
    return(ok)
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   pooled.sd(x, block)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

pooled.sd <- function(x, block) {
    if (!is.numeric(x) && !is(x, "vector") && !is(x, "matrix") && !is(x, "data.frame"))
        stop("The object should be of class numeric, vector, matrix or data.frame! \n")
    if (is(x, "data.frame")) {
        x <- as.matrix(x)
        mode(x) <- "numeric"
    }
    if (any(table(block) == 1))
        stop("All group sample sizes must be greater than 1! \n")
    if (length(table(block)) == 1) {
        if (!is.matrix(x)) {
            return(sd(x, na.rm=TRUE))
        } else {
            return(apply(x, 2, sd, na.rm=TRUE))
        }
    } else {
        block <- as.factor(block)
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        psd <- function(x, tab, n, ng) { sqrt(sum((tab-1)*x^2)/(n - ng)) }
        if (!is.matrix(x)) {
            n <- length(x)
            sdev.gp <- numeric(ng)
            for (g in 1:ng) {
                sdev.gp[g] <- sd(x[def[[g]]], na.rm=TRUE)
            }
            return(psd(x=sdev.gp, tab=tab, n=n, ng=ng))
        } else {
            n <- nrow(x)
            p <- ncol(x)
            sdev.gp <- matrix(data=NA, nrow=ng, ncol=p)
            for (g in 1:ng) {
                sdev.gp[g,] <- apply(x[def[[g]], , drop=FALSE], 2, sd, na.rm=TRUE)
            }
            return(apply(sdev.gp, 2, function(x) psd(x, tab, n, ng)))
        }
    }
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   pooled.mean(x, block)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

pooled.mean <- function(x, block) {
    if (!is.numeric(x) && !is(x, "vector") && !is(x, "matrix") && !is(x, "data.frame"))
        stop("The object should be of class numeric, vector, matrix or data.frame! \n")
    if (is(x, "data.frame")) {
        x <- as.matrix(x)
        mode(x) <- "numeric"
    }
    if (length(table(block)) == 1) {
        if (!is.matrix(x)) {
            return(mean(x, na.rm=TRUE))
        } else {
            return(apply(x, 2, mean, na.rm=TRUE))
        }
    } else {
        block <- as.factor(block)
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        pmean <- function(x, tab, n) { sum(tab*x)/n }
        if (!is.matrix(x)) {
            n <- length(x)
            mean.gp <- numeric(ng)
            for (g in 1:ng) {
                mean.gp[g] <- mean(x[def[[g]]], na.rm=TRUE)
            }
            return(pmean(x=mean.gp, tab=tab, n=n))
        } else {
            n <- nrow(x)
            p <- ncol(x)
            mean.gp <- matrix(data=NA, nrow=ng, ncol=p)
            for (g in 1:ng) {
                mean.gp[g,] <- apply(x[def[[g]], , drop=FALSE], 2, mean, na.rm=TRUE)
            }
            return(apply(mean.gp, 2, function(x) pmean(x, tab, n)))
        }
    }
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   MeanVarReg(data, nc.min, nc.max, probs, B, parallel, conf, verbose)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

MeanVarReg <- function(data, nc.min, nc.max, probs, B, parallel, conf, verbose) {
    data <- t(data)                             # transpose matrix for compatibility with clustering functions
    p <- nrow(data)                             # variables are to be clustered, so they are now by rows !
    n <- ncol(data)
    aver <- apply(data, 1, mean)
    sdev <- apply(data, 1, sd)
    mu.std <- matrix(data=NA, nrow=nc.max - nc.min + 1, ncol=p)
    sd.std <- matrix(data=NA, nrow=nc.max - nc.min + 1, ncol=p)
    gap <- rep(NA, nc.max - nc.min + 1)
    sde <- rep(NA, nc.max - nc.min + 1)
    for (k in nc.min:nc.max) {
        clus <- km.clustering(data=cbind(aver, sdev), k=k, ns=100, maxiter=1000)
        vec.av <- numeric(p)
        vec.sd <- numeric(p)
        for (i in 1:k) {
            wi <- which(clus$membership == i)
            vec.av[wi] <- clus$centers[i,1]
            vec.sd[wi] <- clus$centers[i,2]
        }
        data.std <- as.matrix((data - (vec.av %*% t(rep(1, n)))) / (vec.sd %*% t(rep(1, n))))
        gp <- sim.dis(data=data.std, k=k, B=B, parallel=parallel, conf=conf, verbose=verbose)
        gap[k] <- gp$gapk
        sde[k] <- gp$sk
        mu.std[k,] <- apply(data.std, 1, mean)
        sd.std[k,] <- apply(data.std, 1, sd)
        if (verbose) cat("number of clusters: ", k, " done\n")
    }
    nc <- which.min(gap)
    nc <- max(which(gap <= gap[nc] + sde[nc]))
    clus <- km.clustering(data=cbind(aver, sdev), k=nc, ns=100, maxiter=1000)
    if (!is.null(probs)) {
        mu.quant <- apply(mu.std, 1, function(x) quantile(x=x, probs=probs))
        sd.quant <- apply(sd.std, 1, function(x) quantile(x=x, probs=probs))
    } else {
        mu.quant <- numeric(0)
        sd.quant <- numeric(0)
    }
    return(list(membership=clus$membership, nc=nc, gap=gap, sde=sde,
                mu.std=mu.std, sd.std=sd.std,
                mu.quant=t(mu.quant), sd.quant=t(sd.quant)))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   km.clustering(data, k, ns, maxiter)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

km.clustering <- function(data, k, ns, maxiter) {
    data <- as.matrix(data)
    cn <- unique(data)
    n <- nrow(data)
    p <- ncol(data)
    nn <- nrow(cn)

    Z <- .C("MVR_km_clustering",  as.double(data),
                                  as.double(cn),

                                  centers=double(k*p),
                                  cl=integer(n),
                                  nc=integer(k),
                                  wss=double(k),
                                  tot.wss=as.double(0.0),
                                  err=as.integer(0),

                                  as.integer(n),
                                  as.integer(nn),
                                  as.integer(p),
                                  as.integer(k),
                                  as.integer(ns),
                                  as.integer(maxiter),

                                  NAOK=FALSE,
                                  DUP=TRUE,
                                  PACKAGE="MVR")

    if (Z$err > 0) {
        warning("K-means clustering did not converge in ", maxiter, " iterations.", call.=FALSE)
    }

    cent <- matrix(Z$centers, k)
    dimnames(cent) <- list(1L:k, dimnames(data)[[2L]])

    memb <- Z$cl
    if(!is.null(rn <- rownames(data)))
        names(memb) <- rn

    totss <- sum(scale(data, scale=FALSE)^2)

    obj = structure(list(cluster = memb, centers = cent, totss = totss,
                         withinss = Z$wss, tot.withinss = Z$tot.wss,
                         betweenss = totss - Z$tot.wss, size = Z$nc),
                    class = "kmeans")

    return(list(lWk=log(Z$tot.wss), centers=cent, membership=memb, obj=obj))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   withinsumsq(n, p, B, k)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

withinsumsq <- function(n, p, B, k) {
    W <- .C("MVR_withinsumsq", as.integer(n),
                               as.integer(p),
                               as.integer(k),
                               as.integer(B),
                               lWk=double(B),
                               nstart=as.integer(10),
                               maxiter=as.integer(1000),
                               err=as.integer(0),

                               NAOK=FALSE,
                               DUP=TRUE,
                               PACKAGE="MVR")

    if (W$err > 0) {
        warning("Some k-means clusterings for simulated data did not converge in 1000 iterations.", call.=FALSE)
    }
    return(W$lWk)
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   sim.dis(data, k, B, parallel, conf, verbose)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

sim.dis <- function(data, k, B, parallel, conf, verbose) {
    if (nrow(data) < ncol(data)) data <- t(data)
    lWk <- km.clustering(data=data, k=k, ns=100, maxiter=1000)$lWk
    n <- nrow(data)
    p <- ncol(data)
    if (!parallel) {
        lWk.mc <- withinsumsq(n=n, p=p, B=B, k=k)
    } else {
        if (conf$type == "SOCK") {
            cl <- makeCluster(spec=conf$names,
                              type=conf$type,
                              homogeneous=conf$homo,
                              outfile=conf$outfile,
                              verbose=conf$verbose)
        } else {
            cl <- makeCluster(spec=conf$cpus,
                              type=conf$type,
                              homogeneous=conf$homo,
                              outfile=conf$outfile,
                              verbose=conf$verbose)
        }
        clusterSetRNGStream(cl=cl)
        lWk.cl <- clusterCall(cl=cl, fun=withinsumsq, n=n, p=p, B=ceiling(B/conf$cpus), k=k)
        stopCluster(cl)
        lWk.mc <- numeric(length=0)
        for (j in 1:conf$cpus) {
            lWk.mc <- c(lWk.mc, lWk.cl[[j]])
        }
    }
    lWk.mc.bar <- sum(lWk.mc)/B
    gapk <- abs(lWk.mc.bar - lWk)
    sdk <- sqrt(sum((lWk.mc - lWk.mc.bar)^2)/B)
    return(list(gapk=gapk, sk=sdk * sqrt(1 + (1/B))))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   merging.cluster(M)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

merging.cluster <- function(M) {
    ord <- order(M[, 1, drop=TRUE], decreasing=FALSE)
    N <- M[ord, , drop=FALSE]
    clus <- N[, 1, drop=TRUE]
    if (ncol(M) >= 2) {
        for (g in 2:ncol(M)) {
            memb <- clus
            umemb <- unique(memb)
            clus <- numeric(0)
            i <- 0
            pma <- pmatch(x=names(memb), table=names(N[,g]))
            x <- N[pma,g]
            for (k in umemb) {
                y <- sort(x[memb == k])
                u <- unique(y)
                l <- length(u)
                z <- y
                for (j in 1:l) z[y == u[j]] <- j+i
                i <- i+l
                clus <- c(clus, z)
            }
        }
    }
    return(clus)
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   mvrt(x, obj, lev, tab, ng, def, nc.min, nc.max, parallel, conf, verbose)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

mvrt <- function(x, obj, lev, tab, ng, def, nc.min, nc.max, parallel, conf, verbose) {
    p <- ncol(x)
    aver <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    sdev <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    for (g in 1:ng) {
        aver[g,] <- apply(x[def[[g]], , drop=FALSE], 2, mean, na.rm=TRUE)
        sdev[g,] <- apply(x[def[[g]], , drop=FALSE], 2, sd, na.rm=TRUE)
    }
    M <- matrix(data=NA, nrow=p, ncol=ng, dimnames=list(colnames(x), lev))
    if (is.null(obj)) {
        for (g in 1:ng) {
            if (verbose) cat("Computing optimal cluster configuration of group ", g, " ... \n", sep="")
            M[,g] <- MeanVarReg(data=x[def[[g]], , drop=FALSE], nc.min=nc.min, nc.max=nc.max, probs=NULL, B=100, parallel=parallel, conf=conf, verbose=verbose)$membership
            if (verbose) cat("group : ", g, " done\n")
        }
    } else {
        if (verbose) cat("Retrieving clustering information from 'mvr' object ... \n")
        for (g in 1:ng) {
            M[,g] <- obj$MVR[[g]]$membership
            if (verbose) cat("group : ", g, " done\n")
        }
    }
    clus <- merging.cluster(M)
    aver.reg <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    sdev.reg <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    for (k in unique(clus)) {
        pma <- pmatch(x=colnames(x), table=names(clus))
        ind <- (clus == k)[pma]
        for (g in 1:ng) {
            aver.reg[g,ind] <- mean(aver[g,ind])
            sdev.reg[g,ind] <- sqrt(mean(sdev[g,ind]^2))
        }
    }
    t.reg <- (aver.reg[1,,drop=TRUE] - aver.reg[2,,drop=TRUE]) / sqrt(sdev.reg[1,,drop=TRUE]^2/tab[1] + sdev.reg[2,,drop=TRUE]^2/tab[2])
    return(t.reg)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   .onAttach (libname, pkgname)
################
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("Type MVR.news() to see new features, changes, and bug fixes\n")
}
##########################################################################################################################################
