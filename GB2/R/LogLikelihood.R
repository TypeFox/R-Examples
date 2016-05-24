

# Full GB2 log-likelihood (personal level)
loglp.gb2 <- function(x, shape1, scale, shape2, shape3, w=rep(1, length(x))){
sw <- sum(w)
logf <- logf.gb2(x,shape1,scale,shape2,shape3)
logl <- sum(w*logf)/sw
return(logl)
}

# Full GB2 log-likelihood (household level)
loglh.gb2 <- function (x, shape1, scale, shape2, shape3, w=rep(1, length(x)), hs = rep(1,length(x))) {
    logf <- logf.gb2(x, shape1, scale, shape2, shape3)
    sw <- sum(w*hs)
    loglh <- sum(w*hs*logf)/sw
    return(loglh)
}

# Score functions for the full GB2 log-likelihood (personal level)

scoresp.gb2 <- function (x, shape1, scale, shape2, shape3, w=rep(1, length(x))){ 
    sw <- sum(w)
    lx <- length(x)
    dlogl <- rep(0,4)
    for (i in 1:lx) {
        dlogf <- dlogf.gb2(x[i], shape1, scale, shape2, shape3)
        dlogl <- dlogl + w[i] * dlogf
     }
    return(dlogl/sw)
}

# Score functions for the full GB2 log-likelihood (household level)
scoresh.gb2 <- function(x, shape1, scale, shape2, shape3, w=rep(1, length(x)), hs=rep(1, length(x))){
    sw <- sum(w*hs)
    lx <- length(x)
    dlogl <- rep(0,4)
    for (i in 1:lx) {
        dlogf <- dlogf.gb2(x[i], shape1, scale, shape2, shape3)
        dlogl <- dlogl + w[i] * hs[i]* dlogf
     }
    return(dlogl/sw)
}

# GB2 Fisher information matrix (I_1)
info.gb2 <- function(shape1, scale, shape2, shape3){
I <- matrix(rep(NA,16),ncol=4)
psipq <- digamma(shape2) - digamma(shape3)
trip <- trigamma(shape2)
triq <- trigamma(shape3)
trippq <- trigamma(shape2+shape3)
tripq <- trip + triq
I[1,1] <- (1 + (shape2*shape3/(1+shape2+shape3))*(tripq+(psipq-(shape2-shape3)/(shape2*shape3))^2 - (shape2^2+shape3^2)/(shape2*shape3)^2))/shape1^2
I[1,2] <- (shape2-shape3-shape2*shape3*psipq)/(scale*(1+shape2+shape3))
I[2,1] <- I[1,2]
I[2,2] <- shape1^2*shape2*shape3/(scale^2*(1+shape2+shape3))
I[2,3] <- shape1*shape3/(scale*(shape2+shape3))
I[3,2] <- I[2,3]
I[2,4] <- -shape1*shape2/(scale*(shape2+shape3))
I[4,2] <- I[2,4]
I[1,3] <- -(shape3*psipq-1)/(shape1*(shape2+shape3))
I[3,1] <- I[1,3]
I[3,3] <- trip-trippq
I[1,4] <- (shape2*psipq+1)/(shape1*(shape2+shape3))
I[4,1] <- I[1,4]
I[3,4] <- -trippq
I[4,3] <- -trippq
I[4,4] <- triq-trippq
return(I)
}

