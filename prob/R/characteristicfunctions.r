#  Characteristic functions
#  Released under GPL 2 or greater
#  Copyright January 2009, G. Jay Kerns


cfbeta <- function(t, shape1, shape2, ncp = 0){
    if (shape1 <=0 || shape2 <=0)
        stop("shape1, shape2 must be positive")
    if (identical(all.equal(ncp, 0), TRUE)){
        # require(fAsianOptions)
        kummerM(1i*t, shape1, shape1 + shape2)
    } else {
        fr <- function(x) cos(t*x)*dbeta(x, shape1, shape2, ncp)
        fi <- function(x) sin(t*x)*dbeta(x, shape1, shape2, ncp)
        Rp <- integrate(fr, lower = 0, upper = 1)$value
        Ip <- integrate(fi, lower = 0, upper = 1)$value
        return( Rp + 1i*Ip )
    }
}


cfbinom <- function(t, size, prob){
    if (size <= 0 )
        stop("size must be positive")
    if (prob < 0 || prob > 1)
        stop("prob must be in [0,1]")
    (prob*exp(1i*t) + (1 - prob))^size
}


cfcauchy = function(t, location = 0, scale = 1){
    if (scale <= 0 )
        stop("scale must be positive")
    exp(1i*location*t - scale*abs(t))  
}

cfchisq <- function(t, df, ncp = 0){
    if (df < 0 || ncp < 0  )
        stop("df and ncp must be nonnegative")
    exp(1i*ncp*t/(1-2i*t))/(1 - 2i*t)^(df/2)
}

cfexp <- function(t, rate = 1){
    cfgamma(t, shape = 1, scale = 1/rate)
}

cff <- function(t, df1, df2, ncp, kmax = 10){
    if (df1 <= 0 || df2 <= 0  )
        stop("df1 and df2 must be positive")
    # require(fAsianOptions)
    if( identical(all.equal(ncp, 0), TRUE) ){
        gamma((df1+df2)/2) / gamma(df2/2) * kummerU(-1i*df2*t/df1, df1/2, 1 - df2/2)
    } else {
        exp(-ncp/2)*sum((ncp/2)^(0:kmax)/factorial(0:kmax)* kummerM(-1i*df2*t/df1, df1/2 + 0:kmax, -df2/2))
    }
}


cfgamma <- function(t, shape, rate = 1, scale = 1/rate){
    if (rate <= 0  || scale <= 0)
        stop("rate must be positive")
    (1-scale*1i*t)^(-shape)
}


cfgeom <- function(t, prob){
    cfnbinom(t, size = 1, prob = prob)
}


cfhyper <- function(t, m, n, k){
    if (m < 0 || n < 0 || k < 0)
        stop("m, n, k must be positive")
    hypergeo::hypergeo(-k, -m, n - k + 1, exp(1i*t))/hypergeo::hypergeo(-k, -m, n - k + 1, 1)
}


cflnorm <- function(t, meanlog = 0, sdlog = 1){
    if (sdlog <= 0)
        stop("sdlog must be positive")
    if (identical(all.equal(t, 0), TRUE)){
        return(1+0i)
    } else {
        t <- t * exp(meanlog)
        Rp1 <- integrate(function(y) exp(-log(y/t)^2/2/sdlog^2) * cos(y)/y, lower = 0, upper = t )$value
        Rp2 <- integrate(function(y) exp(-log(y*t)^2/2/sdlog^2) * cos(1/y)/y, lower = 0, upper = 1/t )$value
        Ip1 <- integrate(function(y) exp(-log(y/t)^2/2/sdlog^2) * sin(y)/y, lower = 0, upper = t )$value
        Ip2 <- integrate(function(y) exp(-log(y*t)^2/2/sdlog^2) * sin(1/y)/y, lower = 0, upper = 1/t )$value
        return((Rp1 + Rp2 + 1i*(Ip1 + Ip2))/(sqrt(2*pi) * sdlog))
    }
}


cflogis <- function(t, location = 0, scale = 1){
    if (scale <= 0)
        stop("scale must be positive")
    ifelse( identical(all.equal(t, 0), TRUE),
            return(1),
            return(exp(1i*location)*pi*scale*t/sinh(pi*scale*t)))
}


cfnbinom <- function(t, size, prob, mu){
    if (size <= 0 )
        stop("size must be positive")
    if (prob <= 0 || prob > 1)
        stop("prob must be in (0,1]")
    if (!missing(mu)) {
        if (!missing(prob)) 
            stop("'prob' and 'mu' both specified")
        prob <- size/(size+mu)
    }
    (prob/(1-(1-prob)*exp(1i*t)))^size
}


cfnorm <- function(t, mean = 0, sd = 1){
    if (sd <= 0)
        stop("sd must be positive")
    exp(1i*mean - (sd*t)^2/2)  
}

cfpois <- function(t, lambda){
    if (lambda <= 0)
        stop("lambda must be positive")
    exp(lambda*(exp(1i*t) - 1))
}


cfsignrank <- function(t, n){
    sum(exp(1i*t*0:((n+1)*n/2)) * dsignrank(0:((n+1)*n/2), n))
}


cft <- function(t, df, ncp){
    if(missing(ncp)) ncp <- 0
    if (df <= 0)
        stop("df must be positive")
    if (identical(all.equal(ncp, 0), TRUE)){
        ifelse(identical(all.equal(t, 0), TRUE), 1+0i, 
            as.complex(besselK(sqrt(df)*abs(t), df/2)*(sqrt(df)*abs(t))^(df/2)/( gamma(df/2) * 2^(df/2 - 1) ))
            )
    } else {
        fr <- function(x) cos(t*x)*dt(x, df, ncp)
        fi <- function(x) sin(t*x)*dt(x, df, ncp)
        Rp <- integrate(fr, lower = -Inf, upper = Inf)$value
        Ip <- integrate(fi, lower = -Inf, upper = Inf)$value
        return(Rp + 1i*Ip)
    }
}


cfunif <- function(t, min = 0, max = 1){
    if (max < min)
        stop("min cannot be greater than max")
    ifelse( identical(all.equal(t, 0), TRUE),
            1+0i,
            (exp(1i*t*max) - exp(1i*t*min))/(1i*t*(max - min)))
}


cfweibull <- function(t, shape, scale = 1){
    if (shape <= 0 || scale <= 0)
        stop("shape and scale must be positive")
    fr <- function(x) cos(t*x)*dweibull(x, shape, scale)
    fi <- function(x) sin(t*x)*dweibull(x, shape, scale)
    Rp <- integrate(fr, lower = 0, upper = Inf)$value
    Ip <- integrate(fi, lower = 0, upper = Inf)$value
    return( Rp + 1i*Ip )
}


cfwilcox <- function(t, m, n){
    sum(exp(1i*t*0:(m*n)) * dwilcox(0:(m*n), m, n))
}




