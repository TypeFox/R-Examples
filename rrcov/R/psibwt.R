setMethod("psi", "PsiBwt", function(obj, x) {
    if(obj@c1 == 0){
        x1 <- x - obj@M
        ivec1 <- (x1 <= 0)
        ivec2 <- (x1 > 0)
    } else {
        x1 <- (x-obj@M)/obj@c1
        ivec1 <- (x1 < 0)
        ivec2 <- (x1 >  1)
    }
    return(ivec1*x+(1-ivec1-ivec2)*x*(1-x1^2)^2)
})

setMethod("wt", "PsiBwt", function(obj, x) {
    if(obj@c1 == 0){
        x1 <- x - obj@M
        ivec1 <- (x1 <= 0)
        ivec2 <- (x1 > 0)
    } else {
        x1 <- (x-obj@M)/obj@c1
        ivec1 <- (x1 < 0)
        ivec2 <- (x1 >  1)
    }
    return(ivec1+(1-ivec1-ivec2)*(1-x1^2)^2)
})

setMethod("vt", "PsiBwt", function(obj, x) {
    return(x*psi(obj,x))
})

setMethod("erho", "PsiBwt", function(obj) {
#   expectation of rho(d) under chi-squared p

    c1<-obj@c1
    M<-obj@M
    p<-obj@p

    if (c1 == 0)
            return(.chiInt(p,2,M)/2 + M^2/2*.chiInt2(p,0,M))

    return(.chiInt(p,2,M)/2
        +(M^2/2+c1*(5*c1+16*M)/30)*.chiInt2(p,0,M+c1)
        +(M^2/2-M^2*(M^4-5*M^2*c1^2+15*c1^4)/(30*c1^4))*(.chiInt(p,0,M+c1)-.chiInt(p,0,M))
        +(1/2+M^4/(2*c1^4)-M^2/c1^2)*(.chiInt(p,2,M+c1)-.chiInt(p,2,M))
        +(4*M/(3*c1^2)-4*M^3/(3*c1^4))*(.chiInt(p,3,M+c1)-.chiInt(p,3,M))
        +(3*M^2/(2*c1^4)-1/(2*c1^2))*(.chiInt(p,4,M+c1)-.chiInt(p,4,M))
        -(4*M/(5*c1^4))*(.chiInt(p,5,M+c1)-.chiInt(p,5,M))
        +(1/(6*c1^4))*(.chiInt(p,6,M+c1)-.chiInt(p,6,M)))
})

setMethod("erhoLim", "PsiBwt", function(obj) {
    p<-obj@p
    c1<-obj@c1
    return(.chiInt(p,2,c1) + c1^2*.chiInt2(p,0,c1))
})

setMethod("erhoLimD", "PsiBwt", function(obj) {
    p<-obj@p
    c1<-obj@c1
    return(.chiIntD(p,2,c1) + 2*c1*.chiInt2(p,0,c1)
           + c1^2*.chiInt2D(p,0,c1))
})

setMethod("arpLim", "PsiBwt", function(obj) {
    p<-obj@p
    r<-obj@r
           
    obj@c1 <- c1 <- 2*p
    iter <- 1
    crit <- 100
    eps <- 1e-5
    while(crit>eps & iter<100)
    {
        c1.old <- c1
        fc <- erhoLim(obj) - c1^2*r
        fcp <- erhoLimD(obj) - 2*c1*r
        obj@c1 <- c1 <- c1 - fc/fcp
        if(c1 < 0)  
            obj@c1 <- c1 <- c1.old/2
        crit <- abs(fc)
        ## print(c(iter,c1.old,crit))
        iter <- iter+1
    }

#    print(c(iter,c1,crit))
    return(c(c1, pchisq(c1^2,p), log10(1-pchisq(c1^2,p))))
})

setMethod("csolve", "PsiBwt", function(obj) {
##   find constants c1 and M that give a specified breakdown r
##   and rejection point alpha
    n<-obj@n
    p<-obj@p
    r<-obj@r
    alpha<-obj@alpha
    if(r > (n-p)/(2*n))
        obj@r <- r <- (n-p)/(2*n)         # maximum achievable breakdown

##    print(c(n,p,r,alpha))

    ##
    ##   if rejection is not achievable, use c1=0 and best rejection
    ##
    limvec <- arpLim(obj)
    if(1-limvec[2] <= alpha) {
        obj@c1 <- c1 <- 0
        obj@M <- M <- sqrt(qchisq(1-alpha, p))
##        print("adjusting alpha")
##        print(c(alpha, M, c1))
    }
    else {
        c1.plus.M <- sqrt(qchisq(1-alpha,p))
        obj@M <- M <- sqrt(p)
        obj@c1 <- c1 <- c1.plus.M - M
        iter <- 1
        crit <- 100
        eps <- 1e-5
        while(crit > eps & iter<100){
            deps <- 1e-4
            M.old <- M
            c1.old <- c1
            er <- erho(obj)
            fc <- er - r*(M^2/2+c1*(5*c1+16*M)/30)
            obj1 <- obj; obj1@c1 <- c1+deps
            obj2 <- obj; obj2@M <- M+deps
            fcc1 <- (erho(obj1) - er)/deps
            fcM  <- (erho(obj2) - er)/deps
##            fcp <- fcM - fcc1 - r*(M-(5*c1+16*M)/30+c1*9/30)
            fcp <- fcM - fcc1 - r*(M-(5*c1+16*M)/30+c1*11/30)   # changed according to CB
            obj@M <- M <- M - fc/fcp
            if(M >= c1.plus.M ){
                obj@M <- M <- (M.old + c1.plus.M)/2
            }
            obj@c1 <- c1 <- c1.plus.M - M
            crit <- abs(fc)
##            print(c(iter, M.old, c1.old, crit))
            iter <- iter+1
        }
    }
    return(obj)
})
