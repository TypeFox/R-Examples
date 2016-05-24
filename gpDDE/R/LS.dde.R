make.delay <- function(){
    ## Replace the original proc$dfdc
    dfdc <- function(coefs, bvals, pars, more){
        ## delay <- more$delay
        devals = as.matrix(bvals$bvals %*% coefs)
        ddevals = as.matrix(bvals$dbvals %*% coefs)
        colnames(devals) = more$names
        colnames(ddevals) = more$names
        names(pars) = more$parnames
        ## Need change:
        ## dfdx.d need to be tested.
        g1 <- make.SSElik()$dfdx(ddevals, more$qpts, devals, pars, more)
        g1.d <- more$delay$dfdx.d(ddevals, more$qpts, devals, pars, more)
        weights <- checkweights(more$weights, more$whichobs, g1)
        weights <- mat(weights)
        g2 = weights * (ddevals - more$fn(more$qpts, devals, pars, more$more))
        ## bvals.d <- rbind( matrix(0, dim(g1.d)[1] - dim(more$more$bvals.d)[1], dim(g1.d)[1]), more$more$bvals.d)
        g = as.vector(as.matrix(t(bvals$bvals) %*% g1 + t(more$more$bvals.d) %*% g1.d + 2 * t(bvals$dbvals) %*% g2))
        return(g)
    }


    ## Replacing make.SSElik()$dfdx() for calculating g1.d in the new proc$dfdc
    ## put into proc$more$delay$dfdx.d
    dfdx.d <- function(data, times, devals, pars, more){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights = checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        dfdx.d = more$dfdx.d(times, devals, pars, more$more)
        g = c()
        for (i in 1:dim(dfdx.d)[3]) {
            g = cbind(g, apply(difs * dfdx.d[, , i], 1, sum))
        }
        return(-2 * g)
    }


    ## "more" is taken from $proc$more
    ## To replace the original proc$d2fdc2() has the result of
    ## $$ \bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\bigl[\frac{df_{l}}{dx_{j}dx_{i}}
    ## + \frac{df_{l}}{dx_{j}} \bigl[\frac{df_{l}}{dx_{i}}\bigr]^{\top} $$
    d2fdc2.DDE <- function(coefs, bvals, pars, more){
        delay <- more$delay
        devals = as.matrix(bvals$bvals %*% coefs)
        ddevals = as.matrix(bvals$dbvals %*% coefs)
        colnames(devals) = more$names
        colnames(ddevals) = more$names
        names(pars) = more$parnames
        ## H1: make.SSElik()$d2fdx2 is the same function as lik$d2fdx2 ?!
        ## But the arguments may be different depending on the
        ## Here more is proc$more
        H1 = make.SSElik()$d2fdx2(ddevals, more$qpts, devals, pars,
        more)
        H2 = more$dfdx(more$qpts, devals, pars, more$more)
        weights <- checkweights(more$weights, more$whichobs, mat(H1[,
                                                                    , 1, drop = TRUE]))
        weights <- mat(weights)
##################################################
        ## New for DDE:
##################################################
        H2.1 <- more$dfdx.d(more$qpts, devals, pars, more$more)
        H1.1 <- delay$d2fdx.d2(ddevals, more$qpts, devals, pars, more)
        H1.2 <- delay$d2fdxdx.d(ddevals, more$qpts, devals, pars, more)
        ## H1.3 <- delay$d2fdx.ddx(ddevals, more$qpts, devals, pars, more)
        bvals.d <- rbind( matrix(0, length(H1[,1,1]) - dim(more$more$bvals.d)[1], length(H1[,1,1])), more$more$bvals.d)
        H = list(len = dim(more$more$bvals)[2])
        for (i in 1:dim(devals)[2]){
            H[[i]] = list(len = dim(devals))
            for (j in 1:dim(devals)[2]){
                H[[i]][[j]] <- t(bvals$bvals) %*% diag(H1[, i, j]) %*% bvals$bvals +
                    t(bvals.d) %*% diag(H1.1[, i, j]) %*% bvals.d +
                        t(bvals$bvals) %*% diag(H1.2[, i, j]) %*% bvals.d +
                            t(bvals.d) %*% diag(H1.2[, j, i]) %*% bvals$bvals -
                                2 * t(bvals$dbvals) %*% diag(H2[, i, j] * weights[, i]) %*% bvals$bvals -
                                    2 * t(bvals$bvals) %*% diag(H2[, j, i] * weights[, j]) %*% bvals$dbvals -
                                        2 * t(bvals$dbvals) %*% diag(H2.1[, i, j] * weights[, i]) %*% bvals.d -
                                            2 * t(bvals.d) %*% diag(H2.1[, j, i] * weights[, j]) %*% bvals$dbvals
            }
            H[[i]][[i]] = H[[i]][[i]] + 2 * t(bvals$dbvals) %*% diag(weights[,
                  i]) %*% bvals$dbvals
        }
        H = blocks2mat(H)
        return(H)
    }

##################################################
    ## Used in proc$d2fdc2:
    ##    H2.1 <- more$dfdx.d
    ##    H1.1 <- delay$d2fdx.d2
    ##    H1.2 <- delay$d2fdxdx.d
    ##    H1.3 <- delay$d2fdx.ddx
##################################################

    d2fdx.d2 <- function(data, times, devals, pars, more ){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        dfdx.d = more$dfdx.d(times, devals, pars, more$more)
        d2fdx.d2 = more$d2fdx.d2(times, devals, pars, more$more)
        difs[is.na(difs)] = 0
        weights = checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        H = array(0, c(dim(devals), dim(devals)[2]))
        for (i in 1:dim(d2fdx.d2)[3]) {
            for (j in 1:dim(d2fdx.d2)[4]) {
                H[, i, j] = apply(-difs * d2fdx.d2[, , i, j] + weights *
                 dfdx.d[, , i] * dfdx.d[, , j], 1, sum)
            }
        }
        return(2 * H)
    }

    d2fdxdx.d <- function(data, times, devals, pars, more){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        dfdx = more$dfdx(times, devals, pars, more$more)
        dfdx.d <- more$dfdx.d(times, devals, pars, more$more)
        d2fdxdx.d <- more$d2fdxdx.d(times, devals, pars, more$more)
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        H = array(0, c(dim(devals), dim(devals)[2]))
        for (i in 1:dim(d2fdxdx.d)[3]) {
            for (j in 1:dim(d2fdxdx.d)[4]) {
                H[, i, j] = apply(-difs * d2fdxdx.d[, , i, j] + weights * (dfdx[, , i] * dfdx.d[, , j]), 1, sum)
            }
        }
        return(2 * H)
    }

    d2fdx.ddx <- function(data, times, devals, pars, more){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        dfdx = more$dfdx(times, devals, pars, more$more)
        dfdx.d <- more$dfdx.d(times, devals, pars, more$more)
        d2fdxdx.d <- more$d2fdxdx.d(times, devals, pars, more$more)
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        H = array(0, c(dim(devals), dim(devals)[2]))
        for (i in 1:dim(d2fdxdx.d)[3]){
            for (j in 1:dim(d2fdxdx.d)[4]){
                H[, i, j] = apply(-difs * d2fdxdx.d[, , i, j] + weights *
                 dfdx.d[, , i] * dfdx[, , j], 1, sum)
            }
        }
        return(2 * H)
    }


    ## To replace the original proc$d2fdcdp
    d2fdcdp.sparse <- function (coefs, bvals, pars, more)
    {
        nParsDelay <- sum(more$more$nbeta)
        delay <- more$delay
        devals = as.matrix(bvals$bvals %*% coefs)
        ddevals = as.matrix(bvals$dbvals %*% coefs)
        bvals.d <- more$more$bvals.d
        bvals.d.list <- more$more$bvals.d.list
        ## dbvals.d <- rbind(matrix(0, nrow = dim(bvals$dbvals)[1] - dim(more$more$dbvals.d)[1], ncol = dim(bvals$dbvals)[2]), more$more$dbvals.d)
        colnames(devals) = more$names
        colnames(ddevals) = more$names
        names(pars) = more$parnames
        H1 <- make.SSElik()$d2fdxdp(ddevals, more$qpts, devals, pars,
                                    more)
        H2 = 2 * more$dfdp(more$qpts, devals, pars, more$more)
        weights <- checkweights(more$weights, more$whichobs, mat(H1[,
                                                                    , 1, drop = TRUE]))
        weights <- mat(weights)
##################################################
        ## New for DDE:
##################################################
        H1.1 <- delay$d2fdx.ddp(ddevals, more$qpts, devals, pars, more)
        H3 <- delay$d2fdxdtau(ddevals, more$qpts, devals, pars, more)
        H3.1 <- delay$d2fdx.ddtau(ddevals, more$qpts, devals, pars, more)
        H3.2 <- delay$d2fdx.ddtau2(ddevals, more$qpts, devals, pars, more)
        H4 <- 2 * more$dfdtau(more$qpts, devals, pars, more)
        H = c()
        for (i in 1:(length(pars))){
            H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H1[,, i] + t(bvals.d) %*% H1.1[,,i] - t(bvals$dbvals) %*% (weights * H2[, , i]))))
        }
        for(i in 1: nParsDelay){
            H <- cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H3[,, i] + t(bvals.d) %*% H3.1[,,i] + t(bvals.d.list[[i]]) %*% H3.2[,,i] - t(bvals$dbvals) %*% (weights * H4[,,i]))))
        }
        return(H)
    }

    ## H1.1 <- delay$d2fdx.ddp() :
    ## -2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}d\theta^{\top}}+2\frac{df_{l}}{dx_{j}}\frac{df_{l}}{d\theta^{\top}}
    d2fdx.ddp <- function (data, times, devals, pars, more)
    {
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        dfdx.d = more$dfdx.d(times, devals, pars, more$more)
        dfdp = more$dfdp(times, devals, pars, more$more)
        d2fdx.ddp = more$d2fdx.ddp(times, devals, pars, more$more)
        H = array(0, c(dim(devals), length(pars)-sum(more$tauIndex)))
        for (i in 1:dim(d2fdx.ddp)[3]) {
            for (j in 1:dim(d2fdx.ddp)[4]) {
                H[, i, j] = apply(-difs * d2fdx.ddp[, , i, j] + weights *
                 dfdx.d[, , i] * dfdp[, , j], 1, sum)
            }
        }
        return(2 * H)
    }


    ## H3 <- delay$d2fdxdtau() :
    ## 2\frac{df_{l}}{dx_{j}}\frac{df_{l}}{d\tau^{\top}}-2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}d\tau^{\top}}
    d2fdxdtau <- function (data, times, devals, pars, more)
    {
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        dfdx = more$dfdx(times, devals, pars, more$more)
        dfdtau <- more$dfdtau(times, devals, pars, more)
        d2fdxdtau <- more$d2fdxdtau(times, devals, pars, more)
        H = array(0, c(dim(devals), sum(more$more$nbeta)))
        for (i in 1:dim(d2fdxdtau)[3]) {
            for (j in 1:dim(d2fdxdtau)[4]) {
                H[, i, j] = apply(-difs * d2fdxdtau[, , i, j] + weights *
                 dfdx[, , i] * dfdtau[, , j], 1, sum)
            }
        }
        return(2 * H)
    }

    ## H3.1 <- delay$d2fdx.ddtau() :
    ## 2\frac{df_{l}}{dx_{j}(-T_{j})}\frac{df_{l}}{d\tau^{\top}}-2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}(-T_{j})d\tau^{\top}}
    d2fdx.ddtau <- function(data, times, devals, pars, more){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        dfdx.d = more$dfdx.d(times, devals, pars, more$more)
        dfdtau <- more$dfdtau(times, devals, pars, more)
        d2fdx.ddtau <- more$d2fdx.ddtau(times, devals, pars, more)
        H = array(0, c(dim(devals), sum(more$more$nbeta)))
        for (i in 1:dim(d2fdx.ddtau)[3]) {
            for (j in 1:dim(d2fdx.ddtau)[4]){
                H[, i, j] = apply(-difs * d2fdx.ddtau[, , i, j] + weights *
                 dfdx.d[, , i] * dfdtau[, , j], 1, sum)
            }
        }
        return(2 * H)
    }

    ## H3.2 <- delay$dfdx.ddtau2
    d2fdx.ddtau2 <- function(data, times, devals, pars, more){
        fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights <- checkweights(more$weights, more$whichobs, difs)
        weights <- mat(weights)
        difs = weights * difs
        dfdx.d <- more$dfdx.d(times, devals, pars, more$more)[,,more$more$ndelay, drop = FALSE]
        H = array(0, c(dim(devals), sum(more$more$nbeta)))
        beta.index <- 0
        delay.index <- 0
        for(i in more$more$ndelay){
            delay.index <- delay.index + 1
            for(j in 1:more$more$nbeta[delay.index]){
                beta.index <- beta.index + 1
                H[, i, beta.index] <- apply(-difs * dfdx.d[, , delay.index], 1, sum)
            }
        }
        return(2 * H)
    }
    return(list(dfdc = dfdc, dfdx.d = dfdx.d, d2fdc2.DDE = d2fdc2.DDE, d2fdx.d2 = d2fdx.d2, d2fdxdx.d = d2fdxdx.d, d2fdx.ddx = d2fdx.ddx, d2fdcdp.sparse = d2fdcdp.sparse, d2fdx.ddp = d2fdx.ddp, d2fdxdtau = d2fdxdtau, d2fdx.ddtau = d2fdx.ddtau, d2fdx.ddtau2 = d2fdx.ddtau2))
}

## H4 <- 2 * more$dfdtau() :
## 2\frac{df_{l}}{d\tau^{\top}}
dfdbeta.sparse <- function(times, y, p, more){
    nbeta <- more$more$nbeta
    ndelay <- more$more$ndelay
    x.d <- more$more$y.d.list
    dfdx.d <- more$dfdx.d(times, y, p, more$more)
    r <- array(0, dim = c(dim(dfdx.d)[1:2], sum(nbeta)))
    beta.index <- 0
    delay.index <- 0
    for(j in ndelay){
        delay.index <- delay.index + 1
        for(k in 1:nbeta[delay.index]){
            beta.index <- beta.index + 1
            for(i in 1:dim(dfdx.d)[2]){
                r[,i,beta.index] <- dfdx.d[,i,j] * x.d[[beta.index]]
            }
        }
    }
    ## dimnames(r) <- c(dimnames(dfdx.d)[1:2], dimnames(x.d)[2])
    return(r)
}

d2fxdbeta.sparse <- function(times, y, p, more){
    nbeta <- more$more$nbeta
    ndelay <- more$more$ndelay
    x.d <- more$more$y.d.list
    d2fdxdx.d <- more$d2fdxdx.d(times, y, p, more$more)
    r <- array(0, dim = c(dim(d2fdxdx.d)[1:3], sum(nbeta)))
    beta.index <- 0
    delay.index <- 0
    for(j in ndelay){
        delay.index <- delay.index + 1
        for(k in 1:nbeta[delay.index]){
            beta.index <- beta.index + 1
            for(i in 1:dim(r)[2]){
                for(l in 1:dim(r)[3]){
                    r[,i,l,beta.index] <- d2fdxdx.d[,i,l,j] * x.d[[beta.index]]
                }
            }
        }
    }
    ## dimnames(r) <- c(dimnames(d2fdxdx.d)[1:3], dimnames(x.d))
    return(r)
}

d2fdx.ddbeta.sparse <- function(times, y, p, more){
    nbeta <- more$more$nbeta
    ndelay <- more$more$ndelay
    x.d <- more$more$y.d.list
    d2fdx.d2 <- more$d2fdx.d2(times, y, p, more$more)
    r <- array(0, dim = c(dim(d2fdx.d2)[1:3] ,sum(nbeta)))
    beta.index <- 0
    delay.index <- 0
    for(j in ndelay){
        delay.index <- delay.index + 1
        for(k in nbeta[delay.index]){
            beta.index <- beta.index + 1
            for(i in 1:dim(r)[2]){
                for(l in 1:dim(r)[3]){
                    r[,i,l,beta.index] <- d2fdx.d2[,i,l,j] * x.d[[beta.index]]
                }
            }
        }
    }
    ## dimnames(r) <- dimnames(d2fdx.d2)
    return(r)
}

