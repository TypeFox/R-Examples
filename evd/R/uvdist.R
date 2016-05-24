
"rfrechet"<-
function(n, loc = 0, scale = 1, shape = 1)
{
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc + scale * rexp(n)^(-1/shape)
}

"rgumbel"<-
function(n, loc = 0, scale = 1)
{
    rgev(n, loc = loc, scale = scale, shape = 0)
}

"rrweibull"<-
function(n, loc = 0, scale = 1, shape = 1)
{
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc - scale * rexp(n)^(1/shape)
}

"rnweibull"<-
function(n, loc = 0, scale = 1, shape = 1)
{
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    loc - scale * rexp(n)^(1/shape)
}

"rgev"<-
function(n, loc = 0, scale = 1, shape = 0)
{
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc - scale * log(rexp(n)))
    else return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}

"rgpd"<-
function(n, loc = 0, scale = 1, shape = 0)
{
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(shape == 0) return(loc + scale*rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1) / shape)
}

"rextreme"<-
function(n, quantfun, ..., distn, mlen = 1, largest = TRUE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(largest)
        quantfun(rbeta(n, mlen, 1), ...)
    else
        quantfun(rbeta(n, 1, mlen), ...) 
}

"rorder"<-
function(n, quantfun, ..., distn,  mlen = 1, j = 1, largest = TRUE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen+1-j
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    quantfun(rbeta(n, mlen+1-j, j), ...)
}

"qfrechet"<-
function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc + scale * (-log(p))^(-1/shape)
}

"qgumbel"<-
function(p, loc = 0, scale = 1, lower.tail = TRUE)
{
    qgev(p, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
}

"qrweibull"<-
function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
}

"qnweibull"<-
function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0 || min(shape) <= 0) stop("invalid arguments")
    if(!lower.tail) p <- 1 - p
    loc - scale * (-log(p))^(1/shape)
}

"qgev"<-
function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(!lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale * log(-log(p)))
    else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}

"qgpd"<-
function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(min(scale) < 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    if(lower.tail) p <- 1 - p
    if(shape == 0) return(loc - scale*log(p))
    else return(loc + scale * (p^(-shape) - 1) / shape)
}

"qextreme"<-
function(p, quantfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
{
    if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
        stop("`p' must contain probabilities in (0,1)")
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(quantfun))
        quantfun <- get(paste("q", distn, sep=""), mode="function")
    if(!lower.tail) p <- 1 - p
    if(largest) 
        quantfun(p^(1/mlen), ...)
    else
        quantfun(1-(1-p)^(1/mlen), ...)
}

"pfrechet"<-
function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmax((q - loc)/scale,0)
    p <- exp(-q^(-shape))
    if(!lower.tail) p <- 1 - p
    p
}

"pgumbel"<-
function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    pgev(q, loc = loc, scale = scale, shape = 0, lower.tail = lower.tail)
}

"prweibull"<-
function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    p
}

"pnweibull"<-
function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    q <- pmin((q - loc)/scale,0)
    p <- exp(-(-q)^shape)
    if(!lower.tail) p <- 1 - p
    p
}

"pgev"<-
function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- (q - loc)/scale
    if(shape == 0) p <- exp(-exp(-q))
    else p <- exp( - pmax(1 + shape * q, 0)^(-1/shape))
    if(!lower.tail) p <- 1 - p
    p
}

"pgpd" <-
function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    q <- pmax(q - loc, 0)/scale
    if(shape == 0) p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1/shape)
    }
    if(!lower.tail) p <- 1 - p
    p
}

"pextreme"<-
function(q, distnfun, ..., distn, mlen = 1, largest = TRUE, lower.tail = TRUE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    if(!largest) distn <- 1-distn
    p <- distn^mlen
    if(largest != lower.tail) p <- 1 - p
    p
}

"porder"<-
function(q, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
         lower.tail = TRUE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    lachooseb <- function(a,b) lgamma(a+1) - lgamma(b+1) - lgamma(a-b+1)
    if(largest) svec <- (mlen+1-j):mlen
    else  svec <- 0:(j-1)
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    distn <- distnfun(q, ...)
    store <- matrix(0,nrow=length(q),ncol=j)
    for(k in 1:j)
        store[,k] <- exp(lachooseb(mlen,svec[k]) + svec[k]*log(distn) +
                         (mlen-svec[k])*log(1-distn))
    p <- apply(store,1,sum)
    if(largest != lower.tail) p <- 1 - p
    p
}

"dfrechet"<-
function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    x <- (x - loc)/scale
    xpos <- x[x>0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x>0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x>0 | is.na(x)]
    d <- numeric(nn)
    d[x>0 | is.na(x)] <- log(shape/scale) - (1+shape) * log(xpos) -
         xpos^(-shape)
    d[x<=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dgumbel"<-
function(x, loc = 0, scale = 1, log = FALSE)
{
    dgev(x, loc = loc, scale = scale, shape = 0, log = log)
}

"drweibull"<-
function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    x <- (x - loc)/scale
    xneg <- x[x<0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x<0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x<0 | is.na(x)]
    d <- numeric(nn)
    d[x<0 | is.na(x)] <- log(shape/scale) + (shape-1) * log(-xneg) -
        (-xneg)^shape
    d[x>=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dnweibull"<-
function(x, loc = 0, scale = 1, shape = 1, log = FALSE)
{
    if(min(scale) <= 0 || min(shape) <= 0) stop("invalid arguments")
    x <- (x - loc)/scale
    xneg <- x[x<0 | is.na(x)]
    nn <- length(x)
    scale <- rep(scale, length.out = nn)[x<0 | is.na(x)]
    shape <- rep(shape, length.out = nn)[x<0 | is.na(x)]
    d <- numeric(nn)
    d[x<0 | is.na(x)] <- log(shape/scale) + (shape-1) * log(-xneg) -
        (-xneg)^shape
    d[x>=0 & !is.na(x)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dgev"<-
function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    x <- (x - loc)/scale
    if(shape == 0)
        d <- log(1/scale) - x - exp(-x) 
    else {
        nn <- length(x)
        xx <- 1 + shape*x
        xxpos <- xx[xx>0 | is.na(xx)]
        scale <- rep(scale, length.out = nn)[xx>0 | is.na(xx)]
        d <- numeric(nn)
        d[xx>0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) -
            (1/shape + 1)*log(xxpos)
        d[xx<=0 & !is.na(xx)] <- -Inf
    }  
    if(!log) d <- exp(d)
    d
}

"dgpd"<-
function(x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
    if(min(scale) <= 0) stop("invalid scale")
    if(length(shape) != 1) stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if(shape == 0) {
      d[index] <- log(1/scale[index]) - d[index]
      d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) *
          log(1 + shape * d[index])
        d[!index] <- -Inf
    }
    if(!log) d <- exp(d)
    d
}

"dextreme"<-
function(x, densfun, distnfun, ..., distn, mlen = 1, largest = TRUE, log = FALSE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    if(!largest) distn <- 1 - distn
    distn <- (mlen-1) * log(distn)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- log(mlen) + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
}

"dorder"<-
function(x, densfun, distnfun, ..., distn, mlen = 1, j = 1, largest = TRUE,
         log = FALSE)
{
    if(!is.numeric(mlen) || length(mlen) != 1 || mlen < 1 ||
       mlen %% 1 != 0) 
        stop("`mlen' must be a non-negative integer")
    if(!is.numeric(j) || length(j) != 1 || j < 1 || j %% 1 != 0) 
        stop("`j' must be a non-negative integer")
    if(j > mlen)
        stop("`j' cannot be greater than `mlen'")
    if(!largest) j <- mlen + 1 - j
    if(missing(densfun))
        densfun <- get(paste("d", distn, sep=""), mode="function")
    if(missing(distnfun))
        distnfun <- get(paste("p", distn, sep=""), mode="function")
    dens <- densfun(x, ..., log = TRUE)
    distn <- distnfun(x, ...)[!is.infinite(dens)]
    distn <- (mlen-j) * log(distn) + (j-1) * log(1-distn)
    comb <- lgamma(mlen+1) - lgamma(j) - lgamma(mlen-j+1)
    d <- numeric(length(x))
    d[!is.infinite(dens)] <- comb + dens[!is.infinite(dens)] + distn
    d[is.infinite(dens)] <- -Inf
    if(!log) d <- exp(d)
    d
}





