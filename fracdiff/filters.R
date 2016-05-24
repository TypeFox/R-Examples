#### Provide R implementations of the "obvious" filters needed :

### Have  3 kind of representations:
###
### 1)  (d, AR(p), MA(q))
### 2)  MA(Inf)
### 3)  AR(Inf)


### 'stats' has  ARMAtoMA()  and
###
### arima.sim() which uses filter(*, ) twice  and [after check + 'n.start'] is
##
##    x <- ts(c(rand.gen(n.start, ...), innov[1:n]), start = 1 - n.start)
##    if (length(model$ma)) x <- filter(x, c(1, model$ma), sides = 1)
##    if (length(model$ar)) x <- filter(x, model$ar, method = "recursive")
##    if (n.start > 0)      x <- x[-(1:n.start)]
##    if (d > 0)            x <- diffinv(x, differences = d)
##    as.ts(x)

if(FALSE) { ### Think about case of "arbitrary d" > 0 :
    d0 <- d %% 1  # fractional part
    m <-  d %/% 1 # integer    part
}


fItoMA <- function(d, lag.max)
{
  ## Purpose: ARFIMA(0,d,0) aka {fractional}I(d)  |---> MA representation
  ## ----------------------------------------------------------------------
  ## Arguments: d:       the (fractional) integration order
  ##            lag.max: maximal 'lag' for MA output -- as in ?ARMAtoMA
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 31 May 2004, 16:59

    if(length(d <- as.numeric(d)) != 1)
        stop(sQuote('d')," must be a number {typically in [0,1]}")
    if(length(lag.max <- as.integer(lag.max)) != 1 || lag.max < 1)
        stop(sQuote('lag.max')," must be positive integer")
    psi <- numeric(1+lag.max)
    psi[1] <- 1 # = psi_0
    for(k in 1:lag.max)
        psi[k+1] <- (k-1+d)/k * psi[k]
    psi[-1] # dropping the leading '1'
}

fItoAR <- function(d, lag.max)
{
    ## Purpose: ARFIMA(0,d,0) aka {fractional}I(d)  |---> AR representation
    ## ----------------------------------------------------------------------
    ## Arguments: d:       the (fractional) integration order
    ##            lag.max: maximal 'lag' for MA output -- as in ?ARMAtoMA
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 May 2004, 17:15

    ## Trick:
    fItoMA( -d, lag.max = lag.max)
}


stopifnot(fItoMA(0, 7) == rep(0,7),
          fItoAR(0, 7) == rep(0,7),
          fItoMA(1, 7) == rep(1,7),
          fItoAR(1, 7) == c(-1, rep(0,6)))

### We need  "inf.series" * "inf.series"
###  (or just (1-z)^d  / theta(z)  ---> series
### for the 'residuals' filter  pi(z) = (1-z)^d phi(z) / theta(z)

## package 'polynom' does the following multiplication:
##-                "*" = if(l1 == 1 || l2 == 1) e1 * e2 else {
##-                    m <- outer(e1, e2)
##-                    as.vector(tapply(m, row(m) + col(m), sum))
pMult0 <- function(p1, p2)
{
    ## Purpose: multiplication of polynomials p1 %*% p2
    #### Purpose: multiplication of polynomials (1, p1) %*% (1, p2)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 May 2004, 19:12

    ##m <- outer(c(1,p1), c(1,p2))
    m <- outer(p1, p2)
    as.vector(tapply(m, row(m) + col(m), sum))
}

pMult1 <- function(p1, p2)
{
    ## Purpose: multiplication of polynomials p1 %*% p2   via  FFT
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 May 2004, 19:50
    d1 <- length(p1)-1
    d2 <- length(p2)-1
    Re(fft(fft(c(p1,0[rep(1, d2)])) *
           fft(c(p2,0[rep(1, d1)])), inv = TRUE)) / (d1+d2+1)
}

p1 <- c(1,2)
p2 <- c(2,-1, 3)
(p3 <- pMult0(p1, p2))
(p4 <- pMult0(p2, p3))
stopifnot(p4 == pMult0(p3, p2))

num.Eq <- function(x,y) identical(TRUE, all.equal(x,y))

stopifnot(num.Eq(p3, pMult1(p1, p2)),
          num.Eq(p3, pMult1(p2, p1)),
          num.Eq(p4, pMult1(p2, p3))
          )
