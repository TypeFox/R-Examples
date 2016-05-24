"genmrcheck" <-
function (y, yhat, thresh = -1, sigma = 1, DYADIC = FALSE, beta = 0.5, 
    method = 1, thr.const = 2,schwelle=1000000000) 
{
thresh <- sqrt(thr.const*log(length(y)))
siglevel <- pnorm(thresh)

tmp <-    .C("genmrcheck", y=as.double(y),as.double(yhat), as.integer(length(y)),as.integer(method),as.integer(DYADIC), as.double(siglevel),as.integer(schwelle),as.double(beta),as.double(sigma),PACKAGE="ftnonpar")

tmp$y

}

"genpmreg" <-
function (y, beta = 0.5, squeezing.factor = 0.5, verbose = FALSE, 
    localsqueezing = TRUE, DYADIC = TRUE, thr.const = 2, extrema.nr = -1, 
    bandwidth = -1,SETTOMEAN=FALSE,method=2,...) 
{
    n <- length(y)
    if(method == 2)
      sigma <- mad((y[-1]-y[-n])/sqrt(2))
    else
      sigma <- 1
    #if (bandwidth < 0) 
    #    firstlambda <- 2^(1*floor(log(length(y), base = 2)) - 1)
    if (bandwidth < 0) 
        firstlambda <- 2^(1*floor(log(length(y), base = 2)) - 1)/8
    else firstlambda <- bandwidth
    lambda <- rep(firstlambda, n - 1)
    currprecision <- firstlambda
    while (1 < 2) {
        tmp <- genstring(y, lambda, beta=beta, method=method)
        y.string <- tmp$y
        if(SETTOMEAN)
          {
          if(method == 1)
            y.string <- settomedian(y.string,y)
          else
          if(method == 2)
            y.string <- settomean(y.string,y)
          else
          if(method == 3)
            y.string <- settomean(y.string,y)
          else
          if(method == 4)
            y.string <- settomean(y.string,y)
          else
            print("SETTOMEAN = T might not be a good idea")
          }
        y.mr <- genmrcheck(y,y.string,DYADIC=DYADIC,sigma=sigma,beta=beta,method=method)

        if (verbose) {
            par(mfrow = c(3, 1))
            plot(y, col = "lightgrey",...)
            lines(y.string, col = "red")
            plot(lambda,ty="l",...)
            plot(y.mr,col="green",ty="l")
            print(c("lambda=", min(lambda)))
            print("Press Enter")
            readline()
        }
        if (bandwidth > 0) 
            break
        if (extrema.nr > 0) {
            if (tmp$kext > extrema.nr) 
                lambda <- lambda + currprecision
            if (currprecision < 0.5) {
                if (tmp$kext <= extrema.nr) 
                  break
            }
            else {
                currprecision <- currprecision/2
                lambda <- lambda - currprecision
            }
        }
        else {
            if (sum(y.mr) < 0.1) 
                break
            if (!localsqueezing) 
                lambda <- lambda * squeezing.factor
            else lambda[y.mr[-1] | y.mr[-n]] <- lambda[y.mr[-1] | 
                y.mr[-n]] * squeezing.factor
        }
    }
    list(y = y.string, lambda = lambda, nmax = tmp$kext)
}

"genstring" <-
function (y, lambda, beta = 0.5, method = 1) 
{
    n <- length(y)
    if (length(lambda) == 1) 
        lambda <- c(rep(lambda, n - 1), 0)
    else if (length(lambda) == n - 1) 
        lambda <- c(lambda, 0)
    if (method == 1) {
        tmp <- sort(y)
        eps <- min(c(1e-05 * (max(y) - min(y)), min(tmp[-1] - 
            tmp[-n])/4))
        if (eps < 1e-36) {
            if (mad(y[-1] - y[-n]) > 0) 
                y <- y + rnorm(n, 0, 0.001 * mad(y[-1] - y[-n]))
            else y <- y + rnorm(n, 0, 1e-12)
            tmp <- sort(y)
            eps <- min(tmp[-1] - tmp[-n])/4
        }
    }
    else eps <- 0
    tmp <- .C("genstring", y = as.double(y), as.integer(length(y)), 
        as.double(lambda), as.double(beta), as.integer(method), 
        as.double(eps), kext = as.integer(0),PACKAGE="ftnonpar")
    list(y = tmp$y, kext = tmp$kext)
}


settomedian <-
function(f,y)
{
n <- length(y)

knotsind <- c(1,(2:n)[f[-1] != f[-n]])
knotsy <- f[knotsind]
knotsn <- length(knotsind)

if(knotsn == 1)
  f <- rep(median(y),n)
else
  {
  f[1:(knotsind[2]-1)] <- median(y[1:(knotsind[2]-1)])
  f[knotsind[knotsn]:n] <- median(y[knotsind[knotsn]:n])
  if(knotsn > 2)
    for(i in 2:(knotsn-1))
      {
      if(((knotsy[i]>knotsy[i-1])&&(knotsy[i]>knotsy[i+1]))||((knotsy[i]<knotsy[i-1])&&(knotsy[i]<knotsy[i+1])))
        f[knotsind[i]:(knotsind[i+1]-1)] <- median(y[knotsind[i]:(knotsind[i+1]-1)])
      }
  }
f
}

"l1pmreg" <-
function (y, beta = 0.5, squeezing.factor = 0.5, verbose = FALSE, 
    localsqueezing = TRUE, DYADIC = TRUE, thr.const = 2, extrema.nr = -1, 
    bandwidth = -1,SETTOMEAN=FALSE,method=1,...) 
{
genpmreg(y,beta,squeezing.factor,verbose,localsqueezing,DYADIC,thr.const,
extrema.nr,bandwidth,SETTOMEAN,method,...)
}
"quantpmreg" <- l1pmreg

"settomean" <-
function(f,y)
{
n <- length(y)

knotsind <- c(1,(2:n)[f[-1] != f[-n]])
knotsy <- f[knotsind]
knotsn <- length(knotsind)

if(knotsn == 1)
  f <- rep(mean(y),n)
else
  {
  f[1:(knotsind[2]-1)] <- mean(y[1:(knotsind[2]-1)])
  f[knotsind[knotsn]:n] <- mean(y[knotsind[knotsn]:n])
  if(knotsn > 2)
    for(i in 2:(knotsn-1))
      {
      if(((knotsy[i]>knotsy[i-1])&&(knotsy[i]>knotsy[i+1]))||((knotsy[i]<knotsy[i-1])&&(knotsy[i]<knotsy[i+1])))
        f[knotsind[i]:(knotsind[i+1]-1)] <- mean(y[knotsind[i]:(knotsind[i+1]-1)])
      }
  }
f
}
