`matrices` <-
    function(x,...) UseMethod("matrices")

###{{{ matrices.lvm

mat.lvm <- function(x,ii=index(x),...) {
    A <- ii$A ## Matrix with fixed parameters and ones where parameters are free
    J <- ii$J ## Manifest variable selection matrix
    M1 <- ii$M1 ## Index of free and _unique_ regression parameters
    P <- ii$P  ## Matrix with fixed variance parameters and ones where parameters are free
    P1 <- ii$P1 ## Index of free and _unique_ regression parameters

    constrain.par <- names(constrain(x))
    parval <- list()


    parBelongsTo <- list(mean=seq_len(ii$npar.mean),
                         reg=seq_len(ii$npar.reg)+ii$npar.mean,
                         cov=seq_len(ii$npar.var)+ii$npar.mean+ii$npar.reg,
                         epar=seq_len(ii$npar.ex)+with(ii,npar.reg+npar.var+npar.mean),
                         cpar=numeric())

    idxA <- which(M1==1)
    pidxA <- parBelongsTo$reg
    if (ii$npar.reg>0) {
        A[idxA] <- pidxA
        for (p in ii$parname) {
            idx <- which((x$par==p))
            newval <- A[idx[1]]
            attributes(newval)$reg.idx <- idx
            attributes(newval)$reg.tidx <- which(t(x$par==p))
            parval[[p]] <- newval
            if (length(idx)>1) {
                idxA <- c(idxA,idx[-1])
                pidxA <- c(pidxA,rep(A[idx[1]],length(idx)-1))
                A[idx] <- A[idx[1]]
            }
        } ## duplicate parameters
    }

    pars.var <- parBelongsTo$cov
    idxdiag <- (seq(ncol(P1))-1)*ncol(P1) + seq(ncol(P1))
    idxP <- idxdiag[which(P1[idxdiag]==1)]
    pidxP <- pars.var[seq_len(length(idxP))]
    P[idxP] <- pidxP

    pars.off.diag <- pars.var
    if (length(pidxP)>0) {
        pars.off.diag <- pars.off.diag[-seq_len(length(pidxP))]
    }

    counter <- 0
    if (length(pars.off.diag)>0 & ncol(P)>1)
        for (i in seq_len(ncol(P1)-1))
            for (j in seq(i+1,nrow(P1))) {
                if (ii$P1[j,i]!=0) {
                    counter <- counter+1
                    pos <- c(j+(i-1)*ncol(P1),
                             i+(j-1)*ncol(P1))
                    P[j,i] <- P[i,j] <- pars.off.diag[counter]
                    idxP <- c(idxP,pos); pidxP <- c(pidxP,P[j,i],P[i,j])
                }
            }

    if (length(ii$covparname)>0)
        for (p in ii$covparname) {
            idx <- which(x$covpar==p)
            isOffDiag <- !(idx[1]%in%idxdiag)
            if (!(p%in%ii$parname)) {
                parval[[p]] <- P[idx[1]]
            }
            attributes(parval[[p]])$cov.idx <- idx
            if (length(idx)>1+isOffDiag) {
                P[idx[-seq(1+isOffDiag)]] <- parval[[p]]
            }
            if (ii$npar.reg>0 && p%in%ii$parname) {
                parBelongsTo$reg <- c(parBelongsTo$reg,p)
                idx.reg <- which(x$par==p)
                P[idx] <- A[idx.reg[1]]
                atr <- attributes(parval[[p]])
                parval[[p]] <- A[idx.reg[1]]
                attributes(parval[[p]]) <- atr
                idxP <- c(idxP,idx)
                pidxP <- c(pidxP,rep(P[idx[1]],length(idx)))
            } else {
                idxP <- c(idxP,idx[-seq(1+isOffDiag)])
                pidxP <- c(pidxP,rep(P[idx[1]],length(idx)-1-isOffDiag))
            }
        } ## duplicate parameters

    idxM <- c()
    pidxM <- seq_len(ii$npar.mean)
    v <- NULL

    named <- sapply(x$mean, function(y) is.character(y) & !is.na(y))
    fixed <- sapply(x$mean, function(y) is.numeric(y) & !is.na(y))
    v <- rep(0,length(x$mean))
    names(v) <- colnames(P)
    if (ii$npar.mean>0) {
        idxM <- which(ii$v1==1)
        v[idxM] <- pidxM
    }
    if (any(fixed))
        v[fixed] <- unlist(x$mean[fixed])

    for (p in ii$mparname) {
        idx <- which(x$mean==p)
        if (!(p%in%c(ii$parname,ii$covparname))) {
            if (length(idx)>1) {
                pidxM <- c(pidxM,rep(v[idx[1]],length(idx)-1))
                idxM <- c(idxM,idx[-1])
            }
            parval[[p]] <- v[idx[1]]
            v[idx] <- parval[[p]]
        }
        attributes(parval[[p]])$m.idx <- idx
        if (p %in% ii$covparname & !(p %in% ii$parname)) {
            parBelongsTo$cov <- c(parBelongsTo$cov,p)
            idx.2 <- which(x$covpar==p)
            v[idx] <- P[idx.2[1]]
            pidxM <- c(pidxM,rep(P[idx.2[1]],length(idx)))
            idxM <- c(idxM,idx)
        }
        if (p %in% ii$parname) {
            parBelongsTo$reg <- c(parBelongsTo$reg,p)
            idx.2 <- which(x$par==p)
            v[idx] <- A[idx.2[1]]
            pidxM <- c(pidxM,rep(A[idx.2[1]],length(idx)))
            idxM <- c(idxM,idx)
        }
    }

    ## Ex-parameters
    idxE <- NULL
    pidxE <- parBelongsTo$epar
    named <- sapply(x$exfix, function(y) is.character(y) & !is.na(y))
    fixed <- sapply(x$exfix, function(y) is.numeric(y) & !is.na(y))
    epar <- rep(0,length(x$exfix))
    names(epar) <- names(x$expar)
    if (!(ii$npar.ex==0)) {
        idxE <- which(ii$e1==1)
        epar[idxE] <- pidxE
    }
    if (any(fixed))
        epar[fixed] <- unlist(x$exfix[fixed])
    for (p in ii$eparname) {
        idx <- which(x$exfix==p)
        if (!(p%in%c(ii$parname,ii$covparname,ii$mparname))) {
            if (length(idx)>1) {
                idxE <- c(idxE,idx[-1])
                pidxE <- c(pidxE,rep(epar[idx[1]],length(idx)-1))
            }
            parval[[p]] <- epar[idx[1]]
        }
        attributes(parval[[p]])$e.idx <- idx

        if (length(idx)>1)
            epar[idx[-1]] <- parval[[p]]
        if (p %in% setdiff(ii$covparname,c(ii$parname,ii$mparname))) {
            parBelongsTo$cov <- c(parBelongsTo$cov,p)
            idx.2 <- which(x$covpar==p)
            epar[idx] <- P[idx.2[1]]
            pidxE <- c(pidxE,rep(P[idx.2[1]],length(idx)))
            idxE <- c(idxE,idx)
        }
        if (p %in% setdiff(ii$parname,ii$mparname)) {
            parBelongsTo$reg <- c(parBelongsTo$reg,p)
            idx.2 <- which(x$par==p)
            epar[idx] <- A[idx.2[1]]
            pidxE <- c(pidxE,rep(A[idx.2[1]],length(idx)))
            idxE <- c(idxE,idx)
        }
        if (p %in% ii$mparname) {
            parBelongsTo$mean <- c(parBelongsTo$mean,p)
            idx.2 <- which(x$mean==p)
            epar[idx] <- v[idx.2[1]]
            pidxE <- c(pidxE,rep(v[idx.2[1]],length(idx)))
            idxE <- c(idxE,idx)
        }
    }
    ee <- cbind(idxE,pidxE); rownames(ee) <- names(x$expar)[ee[,1]]

    ## Constrained...
    constrain.par <- names(constrain(x))
    constrain.idx <- NULL
    if (length(constrain.par)>0) {
        constrain.idx <- list()
        for (p in constrain.par) {
            reg.tidx <- reg.idx <- cov.idx <- m.idx <- e.idx <- NULL
            myc <- constrain(x)[[p]]
            xargs <- manifest(x)[na.omit(match(attributes(myc)$args,manifest(x)))]
            if (length(xargs)>0) {
                parval[xargs] <- 0
            }
            if (p%in%ii$parname.all) {
                reg.idx <- which(x$par==p)
                reg.tidx <- which(t(x$par==p))
            }
            if (p%in%ii$covparname.all) {
                cov.idx <- which(x$covpar==p)
            }
            if (p%in%ii$mparname.all) {
                m.idx <- which(x$mean==p)
            }
            if (p%in%ii$eparname.all) {
                e.idx <- which(x$exfix==p)
            }
            constrain.idx[[p]] <- list(reg.idx=reg.idx,reg.tidx=reg.tidx,cov.idx=cov.idx,m.idx=m.idx,e.idx=e.idx)
        }
    }

    parBelongsTo <- lapply(parBelongsTo,function(x) sort(unique(x)))


    return(list(mean=cbind(idxM,pidxM),
                reg=cbind(idxA,pidxA),
                cov=cbind(idxP,pidxP),
                epar=ee,
                parval=parval,
                constrain.idx=constrain.idx,
                parBelongsTo=parBelongsTo))

}



matrices.lvm <- function(x,pars,meanpar=NULL,epars=NULL,data=NULL,...) {
    ii <- index(x)
    pp <- c(rep(NA,ii$npar.mean),pars)
    v <- NULL
    if (!is.null(meanpar) && length(meanpar)>0) {
        pp[seq(ii$npar.mean)] <- meanpar
        v <- ii$v0; v[v==0] <- 0 ##unlist(x$mean[which(v==0)])
        v[ii$mean[,1]] <- meanpar[ii$mean[,2]]
    }
    A <- ii$A
    A[ii$reg[,1]] <- pp[ii$reg[,2]]
    P <- ii$P
    P[ii$cov[,1]] <- pp[ii$cov[,2]]
    e <- NULL
    if (!is.null(epars)) {
        e[ii$epar[,1]] <- pp[ii$epar[,2]]
    }

    parval <- lapply(ii$parval,function(x) {
        res <- pp[x];
        attributes(res) <- attributes(x);
        res })

    ## Constrained...
    constrain.par <- names(constrain(x))
    constrain.idx <- NULL
    cname <- constrainpar <- c()
    if (length(constrain.par)>0 && is.numeric(c(pars,meanpar))) {
        constrain.idx <- list()
        for (p in constrain.par) {
            cname <- c(cname,p)
            myc <- constrain(x)[[p]]
            xargs <- manifest(x)[na.omit(match(attributes(myc)$args,manifest(x)))]
            if (length(xargs)>0) {
                if (!is.null(data)) {
                    parval[xargs] <- (data)[xargs]
                } else parval[xargs] <- 0
            }
            val <- unlist(c(parval,constrainpar,x$mean)[attributes(myc)$args])
            cpar <- myc(val);
            constrainpar <- c(constrainpar,list(cpar)); names(constrainpar) <- cname
            if (p%in%ii$parname.all) {
                if (!is.null(val))
                    A[ii$constrain.idx[[p]]$reg.idx] <- cpar
            }
            if (p%in%ii$covparname.all) {
                if (!is.null(val))
                    P[ii$constrain.idx[[p]]$cov.idx] <- cpar
            }
            if (p%in%ii$mparname.all) {
                if (!is.null(val))
                    v[ii$constrain.idx[[p]]$m.idx] <- cpar
            }
            if (p%in%ii$eparname.all) {
                if (!is.null(val))
                    e[ii$constrain.idx[[p]]$e.idx] <- cpar
            }
        }
    }

    return(list(A=A, P=P, v=v, e=e, parval=parval,
                constrain.idx=ii$constrain.idx, constrainpar=constrainpar))
}

###}}} matrices.lvm

###{{{ matrices.multigroup

matrices.multigroup <- function(x, p, ...) {
  pp <- modelPar(x,p)
  res <- list()
  for (i in seq_len(x$ngroup))
    res <- c(res, list(matrices2(x$lvm[[i]],pp$p[[i]])))
  return(res)
}

###}}}

matrices2 <- function(x,p,...) {
    m0 <- p[seq_len(index(x)$npar.mean)]
    p0 <- p[with(index(x),seq_len(npar)+npar.mean)]
    e0 <- p[with(index(x),seq_len(npar.ex)+npar.mean+npar)]
    matrices(x,p0,m0,e0,...)
}

###{{{ matrices Obsolete
matrices.lvm <- function(x,pars,meanpar=NULL,epars=NULL,data=NULL,...) {
    ii <- index(x)
    A <- ii$A ## Matrix with fixed parameters and ones where parameters are free
    J <- ii$J ## Manifest variable selection matrix
    M0 <- ii$M0 ## Index of free regression parameters
    M1 <- ii$M1 ## Index of free and _unique_ regression parameters
    P <- ii$P  ## Matrix with fixed variance parameters and ones where parameters are free
    P0 <- ii$P0 ## Index of free variance parameters
    P1 <- ii$P1 ## Index of free and _unique_ regression parameters

    P1.lower <- P1[lower.tri(P1)]
    constrain.par <- names(constrain(x))
    parval <- list()

    if (ii$npar.reg>0) {
        A[which(M1==1)] <- pars[seq_len(ii$npar.reg)]
        for (p in ii$parname) {
            idx <- which((x$par==p))
            newval <- A[idx[1]]
            attributes(newval)$reg.idx <- idx
            attributes(newval)$reg.tidx <- which(t(x$par==p))
            parval[[p]] <- newval
            if (length(idx)>1) {
                A[idx[-1]] <- parval[[p]]
            }
        } ## duplicate parameters
    }

    if (ii$npar.reg==0) {
        pars.var <- pars
    } else {
        pars.var <- pars[-seq_len(ii$npar.reg)]
    }
    
    diag(P)[ii$which.diag] <- pars.var[seq_along(ii$which.diag)]

    pars.off.diag <- pars.var
    pars.off.diag <- pars.off.diag[-seq_along(ii$which.diag)]
    counter <- 0
    if (length(pars.off.diag)>0 & ncol(P)>1)
        for (i in seq_len(ncol(P1)-1))
            for (j in seq(i+1,nrow(P1))) {
                if (ii$P1[j,i]!=0) {
                    counter <- counter+1
                    P[j,i] <- pars.off.diag[counter]
                }
            }

    if (length(ii$covparname)>0)
        for (p in ii$covparname) {
            idx <- which(x$covpar==p)
            if (!(p%in%ii$parname)) {
                parval[[p]] <- P[idx[1]]
            }
            attributes(parval[[p]])$cov.idx <- idx
            if (length(idx)>1) {
                P[idx[-1]] <- parval[[p]]
            }
            if (ii$npar.reg>0 && p%in%ii$parname) {
                idx.reg <- which(x$par==p)
                P[idx] <- A[idx.reg[1]]
                atr <- attributes(parval[[p]])
                parval[[p]] <- A[idx.reg[1]] ###?????
                attributes(parval[[p]]) <- atr
            }
        } ## duplicate parameters
    P[upper.tri(P)] <- t(P)[upper.tri(P)]  ## Symmetrize...


    v <- NULL
    {
        named <- sapply(x$mean, function(y) is.character(y) & !is.na(y))
        fixed <- sapply(x$mean, function(y) is.numeric(y) & !is.na(y))
        v <- rep(0,length(x$mean))
        names(v) <- colnames(P)
        if (!(is.null(meanpar) | ii$npar.mean==0))
            v[ii$v1==1] <- meanpar
        if (any(fixed))
            v[fixed] <- unlist(x$mean[fixed])

        for (p in ii$mparname) {
            idx <- which(x$mean==p)

            if (!(p%in%c(ii$parname,ii$covparname))) {
                parval[[p]] <- v[idx[1]]
            }
            attributes(parval[[p]])$m.idx <- idx

            if (length(idx)>1)
                v[idx[-1]] <- parval[[p]]
            if (p %in% ii$covparname & !(p %in% ii$parname)) {
                idx.2 <- which(x$covpar==p)
                v[idx] <- P[idx.2[1]]
            }
            if (p %in% ii$parname) {
                idx.2 <- which(x$par==p)
                v[idx] <- A[idx.2[1]]
            }
        }
    }


    ## Ex-parameters
    e <- NULL
    {
        named <- sapply(x$exfix, function(y) is.character(y) & !is.na(y))
        fixed <- sapply(x$exfix, function(y) is.numeric(y) & !is.na(y))

        e <- rep(0,length(x$exfix))

        names(e) <- names(x$expar)
        if (!(is.null(epars) | ii$npar.ex==0))
            e[which(ii$e1==1)] <- epars
        if (any(fixed))
            e[fixed] <- unlist(x$exfix[fixed])
        for (p in ii$eparname) {
            idx <- which(x$exfix==p)
            if (!(p%in%c(ii$parname,ii$covparname,ii$mparname))) {
                parval[[p]] <- e[idx[1]]
            }
            attributes(parval[[p]])$e.idx <- idx

            if (length(idx)>1)
                e[idx[-1]] <- parval[[p]]
            if (p %in% setdiff(ii$covparname,c(ii$parname,ii$mparname))) {
                idx.2 <- which(x$covpar==p)
                e[idx] <- P[idx.2[1]]
            }
            if (p %in% setdiff(ii$parname,ii$mparname)) {
                idx.2 <- which(x$par==p)
                e[idx] <- A[idx.2[1]]
            }
            if (p %in% ii$mparname) {
                idx.2 <- which(x$mean==p)
                e[idx] <- v[idx.2[1]]
            }
        }
    }

    ## Constrained...
    constrain.idx <- NULL
    cname <- constrainpar <- c()
    if (length(constrain.par)>0 && is.numeric(c(pars,meanpar,e))) {
        constrain.idx <- list()
        for (p in constrain.par) {
            cname <- c(cname,p)
            reg.tidx <- reg.idx <- cov.idx <- m.idx <- e.idx <- NULL
            myc <- constrain(x)[[p]]
            xargs <- manifest(x)[na.omit(match(attributes(myc)$args,manifest(x)))]
            if (length(xargs)>0) {
                if (!is.null(data)) {
                    parval[xargs] <- (data)[xargs]
                } else parval[xargs] <- 0
            }
            val <- unlist(c(parval,constrainpar,x$mean,e)[attributes(myc)$args])
            cpar <- myc(val);
            constrainpar <- c(constrainpar,list(cpar)); names(constrainpar) <- cname
            if (p%in%ii$parname.all) {
                reg.idx <- which(x$par==p)
                reg.tidx <- which(t(x$par==p))
                if (!is.null(val))
                    A[reg.idx] <- cpar##myc(val)
            }
            if (p%in%ii$covparname.all) {
                cov.idx <- which(x$covpar==p)
                if (!is.null(val))
                    P[cov.idx] <- cpar##myc(val)
            }
            if (p%in%ii$mparname.all) {
                m.idx <- which(x$mean==p)
                if (!is.null(val))
                    v[m.idx] <- cpar##myc(val)
            }
            if (p%in%ii$eparname.all) {
                e.idx <- which(x$exfix==p)
                if (!is.null(val))
                    e[e.idx] <- cpar##myc(val)
            }
            constrain.idx[[p]] <- list(reg.idx=reg.idx,reg.tidx=reg.tidx,cov.idx=cov.idx,m.idx=m.idx,e.idx=e.idx)
        }
    }

    if (x$index$sparse & !is.character(class(pars)[1])) {
        A <- as(A,"sparseMatrix")
        P <- as(P,"sparseMatrix")
        v <- as(v,"sparseMatrix")
    }
    return(list(A=A, P=P, v=v, e=e, parval=parval, constrain.idx=constrain.idx, constrainpar=constrainpar))
}
###}}} matrices Obsolete
