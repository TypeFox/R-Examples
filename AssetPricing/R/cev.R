cev <- function(x,t,v,type,maximize=FALSE) {
#
# Conditional expected value (given an arrival at
# the time in question).  Note that objects dS, gpr, alpha
# and epsilon are assigned in the environment of cev.
#
if(!(type%in%c("sip","dip")))
    stop(paste("Type",type,"not recognized.\n"))
qmax <- length(v)
jmax <- length(gpr)
R    <- numeric(qmax)

if(maximize) {
# Since we are maximizing we are using either xsolve.disc() or
# xsolve.pwl().  If we are using xsolve.disc() then x is a vector
# whose entries are the possible discrete prices.  If we are using
# xsolve.pwl() then x is a list, created by getPossPrices(), with
# one entry for each value of q.
        E      <- parent.env(environment())
        xback  <- E$xback
        usexb  <- !is.null(xback)
        lopt   <- if(type=="sip") qmax else jmax*(qmax+0.5-jmax/2)
        xopt   <- numeric(lopt)
        for(q in 1:qmax) {
            vq   <- v[q]
            jtop <- min(q,jmax)
            rv   <- (c(v[q:1],0)[-1])[1:jtop]
            K    <- gpr[1:jtop]
            if(jtop < jmax) K[jtop] <- K[jtop] + alpha*(1-sum(K))
            xc <- if(is.list(x)) x[[q]] else x
            if(type=="sip") {
                if(usexb) {
                    xb <- xback[q]
                    xc <- sort(unique(c(xc,xb)))
                }
                nx  <- length(xc)
                tmp <- numeric(nx)
                S   <- getS(dS,xc,t,1:jtop)
                for(i in 1:nx) {
                    tmp[i] <- xc[i]*sum((1:jtop)*S[i,1:jtop]*K) +
                                    sum((rv - vq)*S[i,1:jtop]*K)
                }
# New stuff.  If there is "no discernable improvement" by changing
# to a new price, i.e. if the maximum of tmp is no more than its
# value at "the previous price" minus epsilon, then set the new
# price equal to the previous value.  If there is "discernable
# improvement" then look at those values of price where the value
# of tmp is greater than its maximum value minus epsilon and pick
# the one which is closest to the previous price.  This prevents
# there being too many jumps, and where there are jumps, prevents
# them from being unnecessarily large.  Before, we just took the
# first of the possible maxima of tmp.  If a number of values of
# tmp differ only by "numerical noise" then the first maximum can
# occur at a price substantially different from the previous price
# when there is no meaningful improvement from the change.
                if(usexb) {
                    mt <- max(tmp)
                    ib <- which(xc==xb)
                    if(mt - tmp[ib] <= epsilon) im <- ib else {
                        ic <- which(tmp > mt - epsilon)
                        im <- ic[which.min(abs(xc[ic]-xb))]
                    }
                } else im <- which.max(tmp)[1]
# End new stuff.
                R[q]    <- tmp[im] + vq
                xopt[q] <- xc[im]
                } else { # type == "dip"
# Here xc may be a vector of candidate prices --- discrete pricing ---
# or a list whose j-th entry is the vector of candidate prices for
# groups of size j.
                    iqj <- qj2i(q,1:jtop,qmax)
                    if(usexb) xb  <- xback[iqj]
                    Rtmp <- numeric(jtop)
                    xtmp <- numeric(jtop)
                    for(j in 1:jtop) {
                        xcj <- if(is.list(xc)) xc[[j]] else xc
                        if(usexb) xcj <- sort(unique(c(xcj,xb[j])))
                        nxj <- length(xcj)
                        S   <- getS(dS,xcj,t,j)
                        for(i in 1:nxj) {
                            tmp[i] <- (j*xcj[i] + rv[j] - vq)*S[i,1]*K[j]
                        }
# New stuff for the "dip" case.
                        if(usexb) {
                            mt <- max(tmp)
                            ib <- which(xcj==xb[j])
                            if(mt - tmp[ib] <= epsilon) im <- ib else {
                                ic <- which(tmp > mt - epsilon)
                                im <- ic[which.min(abs(xcj[ic]-xb))]
                            }
                        } else im  <- which.max(tmp)[1]
# End new stuff.
                        Rtmp[j] <- tmp[im]
                        xtmp[j] <- xc[im]
                    }
                R[q] <- sum(Rtmp) + vq
                iv   <- qj2i(q,1:jtop,qmax)
                xopt[iv] <- xtmp
            }
        }
    attr(R,"xopt") <- xopt
} else {
S   <- getS(dS,x,t,1:jmax)
    for(q in 1:qmax) {
        vq   <- v[q]
        jtop <- min(q,jmax)
        rv   <- (c(v[q:1],0)[-1])[1:jtop]
        K    <- gpr[1:jtop]
        if(jtop < jmax) K[jtop] <- K[jtop] + alpha*(1-sum(K))
        if(type=="sip") {
            R[q] <- x[q]*sum((1:jtop)*S[q,1:jtop]*K) +
                                 sum((rv - vq)*S[q,1:jtop]*K) + vq
        } else {
            i    <- qj2i(q,1:jtop,qmax)
            R[q] <- sum(((1:jtop)*x[i] + rv - vq)*
                                        S[cbind(i,1:jtop)]*K) + vq
        }
    }
}
R
}
