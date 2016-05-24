# Pricing of different certificate types : Bonus certificate, Outperformance etc.
#
# Stefan Wilhelm

#library(fOptions)
#library(fExoticOptions)

# Pricing of a Warrant (Options)
#
# @params S underlying price
# @params X strike price of the embedded option
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
Warrant<-function(type, S, X, Time, r, r_d, sigma, ratio=1)
{
  if (Time == 0)
  {
    if (type == "c")
    {
      pmax(S-X, 0) * ratio
    } 
    else
    {
      pmax(X-S, 0) * ratio
    }
  }
  else
  {
    option <- GBSOption(TypeFlag=ifelse(type == "c","c","p"), S, X, Time, r, b=r-r_d, sigma) 
    price1 <- pmax(attr(option,"price"),0)
    price1 * ratio
  }
}


# Pricing of a Straddle
#
# @params S underlying price
# @params X strike price of the embedded option
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
Straddle<-function(S, X, Time, r, r_d, sigma, ratio=1)
{
  if (Time == 0)
  {
    (pmax(S-X, 0) + pmax(X-S, 0)) * ratio
  }
  else
  {
    call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma) 
    put <- GBSOption(TypeFlag="p", S, X, Time, r, b=r-r_d, sigma) 
    price1 <- pmax(attr(call,"price"),0)
    price2 <- pmax(attr(put,"price"),0)
    (price1 + price2) * ratio
  }
}

# Pricing of TurboCertificate
#
# Duplication:
# (1) Barrier Options
#
# @params S underlying price
# @params X strike price of the embedded option
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
TurboCertificate<-function(type, S, X, B, Time, r, r_d, sigma, ratio=1)
{
  if (Time == 0)
  {
    if (type == "c")
    {
      ifelse (S>=B , pmax(S-X, 0), 0) * ratio
    }
    else
    {
      ifelse (S<=B , pmax(X-S, 0), 0) * ratio
    }
  }
  else
  {
    barrieroption = StandardBarrierOption(TypeFlag=ifelse(type == "c","cdo","puo"), S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price1 <- pmax(attr(barrieroption,"price"),0)
    price1 * ratio
  }
}

# Pricing of Strangle
#
# Duplication:
# (1) Long-Call mit Strike X2
# (2) Long-Put mit Strike X1
#
Strangle <- function (S, X1, X2, Time, r, r_d, sigma, ratio = 1) 
{
    if (Time == 0) {
        (pmax(S - X2, 0) + pmax(X1 - S, 0)) * ratio
    }
    else {
        call <- GBSOption(TypeFlag = "c", S, X2, Time, r, b = r - 
            r_d, sigma)
        put <- GBSOption(TypeFlag = "p", S, X1, Time, r, b = r - 
            r_d, sigma)
        price1 <- pmax(attr(call, "price"), 0)
        price2 <- pmax(attr(put, "price"), 0)
        (price1 + price2) * ratio
    }
}

# Pricing of TurboStrangle
#
# Duplication:
# (1) Down-And-Out-Call
# (2) Up-And-Out-Put
#
# @params S underlying price
# @params X.Call strike price of the embedded barrier call option
# @params B.Call barrier of the embedded barrier call option
# @params X.Put strike price of the embedded barrier put option
# @params B.Put barrier of the embedded barrier put option
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
TurboStrangle<-function(S, X.Call, B.Call, X.Put, B.Put, Time, r, r_d, sigma, ratio=1)
{
  if (Time == 0)
  {
    (ifelse (S>=B.Call , pmax(S-X.Call, 0), 0) + ifelse (S<=B.Put, pmax(X.Put-S, 0), 0)) * ratio
  }
  else
  {
    turbo_call <- StandardBarrierOption(TypeFlag="cdo", S, X=X.Call, H=B.Call, K=0, Time, r, b=r-r_d, sigma) 
    price1 <- pmax(attr(turbo_call,"price"),0)
    #price1=0
    
    turbo_put <- StandardBarrierOption(TypeFlag="puo", S, X=X.Put, H=B.Put, K=0, Time, r, b=r-r_d, sigma) 
    price2 = pmax(attr(turbo_put,"price"),0)
    #price2=0
    (price1 + price2) * ratio
  }
}

#######################################################
#
# Berechnung der Griechen
#
#######################################################

# Berechne das Delta für eine Preisfunktion über den Differenzenquotienten
.Delta<-function(FUN, S, dS=0.001, ...)
{
  p1<-FUN(S=S, ...)
  p2<-FUN(S=S+dS, ...)
  delta<-(p2-p1)/dS
  delta 
}

# Berechne das Vega für eine Preisfunktion über den Differenzenquotienten
.Vega<-function(FUN, sigma, dv=0.001, ...)
{
  p1<-FUN(sigma=sigma, ...)
  p2<-FUN(sigma=sigma+dv, ...)
  vega<-(p2-p1)/dv
  vega 
}

# Berechne das Rho für eine Preisfunktion über den Differenzenquotienten
.Rho<-function(FUN, r, dr=0.001, ...)
{
  p1<-FUN(r=r, ...)
  p2<-FUN(r=r+dr, ...)
  rho<-(p2-p1)/dr
  rho 
}

# Berechne das Theta für eine Preisfunktion über den Differenzenquotienten
.Theta<-function(FUN, Time, dT=0.001, ...)
{
  p1<-FUN(Time=Time, ...)
  p2<-FUN(Time=Time+dT, ...)
  theta<-(p2-p1)/dT
  theta 
}


# Pricing of DiscountCertificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Short Call
#
# @params S underlying price
# @params X strike price of the embedded option (cap)
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
DiscountCertificate<-function(S, X, Time, r, r_d, sigma, ratio=1)
{
  if (Time == 0)
  {
    pmax(pmin(S,X),0) * ratio
  }
  else
  {
    # 1. Zero-Strike Call
    zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
    price1 <- attr(zero_strike_call,"price")
    
    # 2. European Short-Call with Strike X (=Cap)
    short_call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma) 
    price2 <- attr(short_call,"price")
    
    (price1 - price2) * ratio
  }
}

# Eindimensionales Newton-Verfahren mit Differenzenquotienten statt Ableitung
## Intervall wird hier fuer Demo nur zum Zeichnen verwendet
#
# @param f Funktion f(x)
# @param x Startwert
newton <- function(f, x, interval, tol=sqrt(.Machine$double.eps), doPlot=FALSE)
{
  a <- min(interval)
  b <- max(interval)
  
  if (doPlot)
  {
    plot(f, a, b, lwd=2)
    abline(h=0, lty=2)
  }
  n <- 0
  fx <- Inf
  while(abs(fx)>tol)
  {
   fx <- f(x)
   # Differenzenquotient
   gx <- (fx-f(x-0.001))/0.001
   		
   if (doPlot)
   {		
    abline(a=fx-x*gx, b=gx, col=3)  # Tangente
   }
      
   x <- x - fx/gx
   n <- n+1
   if (x < 0) break
  }
  list(root=x, f.root=fx, iter=n)
}

# imply volatility for given price function f and a given price "price" using one-dimensional Newton-Raphson
#
implyVolatility <- function(price, f, interval=c(0, 1), sigma=NULL, doPlot=FALSE, ...)
{
  if (is.null(sigma) || is.infinite(sigma)) sigma = 0.2
  newton(f=function(sigma){f(sigma=sigma, ...) - price}, x=sigma, interval=interval, doPlot=doPlot)
}

# Pricing of DiscountPlusCertificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Short Call
# (3) Down-And-Out-Put
#
# @params S underlying price
# @params X strike price of the embedded option (cap)
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
# @params barrierActive
# @param barrierHit
DiscountPlusCertificate<-function(S, X, B, Time, r, r_d, sigma, ratio=1, barrierActive=TRUE, barrierHit=FALSE)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. European Short-Call with Strike X (=Cap)
  short_call     = GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma) 
  price2 <- attr(short_call,"price")
  price2 <- ifelse(is.na(price2), 0, price2)
  
  if (barrierActive && !barrierHit)
  {
    # 3. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
    down_out_put <- StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price3 <- pmax(attr(down_out_put,"price"),0)
	price3=ifelse(is.na(price3), 0, price3)
  }
  else
  {
    price3 <- 0
  }
  
  (price1 - price2 + price3) * ratio
}

# Pricing of ReverseDiscountCertificate (Reverse Discount)
#
# Duplication:
# (1) Short Position in Aktie
# (2) European Short-Put
#
# @param S   underlying price
# @param S0  Startkurs
# @param X   strike price of the embedded option (cap)
# @param Time time to maturity in years
# @param r interest rate p.a. as 0.02 = 2%
# @param r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @param ratio
ReverseDiscountCertificate<-function(S, S0, X, Time, r, r_d, sigma, ratio=1)
{
  # 1. Short-Position in Aktie
  #price1=2*S0 - S
  put <- GBSOption(TypeFlag="p", S, X=2*S0, Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(put,"price"),0)
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. European Short-Put with Strike X (=Cap)
  short_put <- GBSOption(TypeFlag="p", S, X, Time, r, b=r-r_d, sigma) 
  price2 <- attr(short_put,"price")
  price2 <- ifelse(is.na(price2), 0, price2)
  
  (price1 - price2) * ratio
}


# Pricing of ReverseDiscountPlusCertificate (Reverse Protect Discount Plus Pro)
#
# Duplication:
# (1) Short Position in Aktie
# (2) European Short-Put
# (3) Up-And-Out-Call with Strike X and Barrier B
#
# @param S underlying price
# @param S0 Startkurs
# @param X strike price of the embedded option (cap)
# @param B barrier
# @param Time time to maturity in years
# @param r interest rate p.a. as 0.02 = 2%
# @param r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @param ratio
# @param barrierActive
ReverseDiscountPlusCertificate<-function(S, S0, X, B, Time, r, r_d, sigma, ratio=1, barrierActive=TRUE)
{
  # 1. Short-Position in Aktie
  #price1=2*S0 - S
  put <- GBSOption(TypeFlag="p", S, X=2*S0, Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(put,"price"),0)
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. European Short-Put with Strike X (=Cap)
  short_put <- GBSOption(TypeFlag="p", S, X, Time, r, b=r-r_d, sigma) 
  price2 <- attr(short_put,"price")
  price2 <- ifelse(is.na(price2), 0, price2)
  
  if (barrierActive)
  {
    # 3. Up-And-Out-Call with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
    up_out_call <- StandardBarrierOption(TypeFlag="cuo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price3 <- pmax(attr(up_out_call,"price"),0)
    price3 <- ifelse(is.na(price3), 0, price3)
  }
  else
  {
    price3 <- 0
  }
  
  (price1 - price2 + price3) * ratio
}

# Pricing of DiscountCall
#
# Duplication:
# (1) Long-Call
# (2) Short-Call
#
# @params S underlying price
# @params X strike price of the embedded long Call option (strike)
# @params Cap strike price of the embedded Short Call option (Cap)
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
DiscountCall<-function(S, X, Cap, Time, r, r_d, sigma, ratio=1)
{
  # 1. European Long-Call with Strike X (=Cap)
  if (Time == 0)
  {
    price1 <- pmax(S-X,0)
  }
  else
  {
    long_call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma) 
    price1 <- attr(long_call,"price")
    price1 <- ifelse(is.na(price1), 0, price1)
  }
  
  # 2. European Short-Call with Strike Cap (=Cap)
  if (Time == 0)
  {
    price2 <- pmax(S-Cap,0)
  }
  else
  {
    short_call <- GBSOption(TypeFlag="c", S, Cap, Time, r, b=r-r_d, sigma) 
    price2 <- attr(short_call,"price")
    price2 <- ifelse(is.na(price2), 0, price2)
  }
  
  (price1 - price2) * ratio
}


# Pricing of DiscountPut
#
# Duplication:
# (1) Long-Put mit Strike = Cap
# (2) Short-Put mit Strike = X (X < Cap)
#
# @params S underlying price
# @params X strike price of the embedded long Put option (strike)
# @params Cap strike price of the embedded Short Put option (Cap)
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params ratio
DiscountPut<-function(S, X, Cap, Time, r, r_d, sigma, ratio=1)
{
  # 1. European Long-Put with Strike =Cap
  if (Time == 0)
  {
    price1 = pmax(Cap-S,0)
  }
  else
  {
    
    long_put <- GBSOption(TypeFlag="p", S, Cap, Time, r, b=r-r_d, sigma) 
    price1 <- pmax(attr(long_put,"price"),0)
  }
  
  # 2. European Short-Put with Strike =X
  if (Time == 0)
  {
    price2 <- pmax(X-S,0)
  }
  else
  {
    short_put <- GBSOption(TypeFlag="p", S, X, Time, r, b=r-r_d, sigma) 
    price2 <- pmax(attr(short_put,"price",0))
  }
  
  (price1 - price2) * ratio
}



# Pricing of BonusCertificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Down-And-Out-Put
#
# @params S underlying price
# @params X strike price (bonus level)
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
BonusCertificate<-function(S, X, B, Time, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  
  # 2. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  if (barrierHit)
  {
    price2 <- 0
  }
  else
  {
    down_out_put <- StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price2 <- pmax(attr(down_out_put,"price"),0)
	  price2 <- ifelse(is.na(price2), 0, price2)
  }
  #price2 = ifelse(barrierHit, 0, pmax(attr(StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma),"price"),0)
  
  (price1 + price2) * ratio
}

# Bonus-Pro-Zertifikat: Barriere ist erst ab Zeitpunkt t1<=Time aktiv.
#
# Modellierung durch eine Partial-Time-Barrier-Option (PTSingleAssetBarrierOption) Typ "pdoB2"
#
# Literatur: 
# 
# - Haug (2007), S.160 ff.
# - Heynen & Kat (1994)
#
# @param TypeFlag
# There are two types of "B" options: 
# "B1" is defined such that only a barrier hit or crossed causes the option to be knocked out, and a 
# "B2" is defined such that a down-and-out-put is knocked out as soon as the underlying price is below the barrier.
#
# TypeFlag = "poB1" : 
# The barrier of the down-and-out-put is only monitored in [time1, Time] with 0 <= time1 <= Time (partial-time monitoring) 
# instead of [0, Time]. Ceteris paribus, this means a reduced risk of knock-out.
# For time1 = 0 (full-time monitoring), the value of a Bonus Pro equals the value of a standard Bonus certificate.
# For time1 = Time (no barrier to be monitored), the value of the type "poB1" Bonus Pro duplicates a Protective Put strategy (except for the dividend payments).
#
# TypeFlag = "pdoB2":
#
# @param time1 start time of barrier monitoring
#
BonusProCertificate <- function(TypeFlag=c("poB1","pdoB2"), S, X, B, Time, time1=0, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  TypeFlag = match.arg(TypeFlag)
  
  if (time1 == 0) {
    p <- BonusCertificate(S=S, X=X, B=B, Time=Time, r=r, r_d=r_d, sigma=sigma, ratio=ratio, barrierHit=barrierHit)
    return(p)
  }
  
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  
  # 2. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  price2 <- ifelse (barrierHit | (time1 == 0 & S <= B), 0,  pmax(PTSingleAssetBarrierOption(TypeFlag, S, X, H=B, time1=time1, Time2=Time, r, b=r-r_d, sigma)@price,0))
  
  
  (price1 + price2) * ratio
}

# Pricing of Capped Bonus Certificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Down-And-Out-Put
# (3) Short Call
#
# @params S underlying price
# @params X strike price (bonus level)
# @params B barrier
# @params Cap cap
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params ratio
CappedBonusCertificate<-function(S, X, B, Cap, Time, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  
  # 2. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  price2 <- ifelse (barrierHit | (S < B), 0, pmax(StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma)@price, 0))
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 3. European Short-Call with Strike X (=Cap)
  short_call <- GBSOption(TypeFlag="c", S, X=Cap, Time, r, b=r-r_d, sigma) 
  price3 <- attr(short_call,"price")
  price3 <- ifelse(is.na(price3), 0, price3)
  
  (price1 + price2 - price3) * ratio
}

# Pricing of Reverse Bonus Certificate
#
# Duplication:
# (1) Short Position in Aktie bezogen auf S0
# (2) Long Up-And-Out-Call mit Strike X und Barrier B
#
# @params S underlying price
# @params S0 underlying start price
# @params X strike price (bonus level)
# @params B barrier
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
ReverseBonusCertificate<-function(S, S0, X, B, Time, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  # 1. Short-Position in Aktie : entspricht einem Long-Put mit Strike 2*S0
  ##price1= pmax(2*S0*exp(-T*(r-r_d)) - S,0)
  put <- GBSOption(TypeFlag="p", S, X=2*S0, Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(put,"price"),0)
  
  if (barrierHit)
  {
    price2 = 0
  }
  else
  {
    # 2. Up-And-Out-Call with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
    up_out_call <- StandardBarrierOption(TypeFlag="cuo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price2 <- pmax(attr(up_out_call,"price"),0)
    price2 <- ifelse(is.na(price2), 0, price2)
  }
  
  (price1 + price2) * ratio
}

# Pricing of Capped Reverse Bonus Certificate
#
# Duplication:
# (1) Short Position in Aktie bezogen auf S0
# (2) Long Up-And-Out-Call mit Strike X und Barrier B
# (3) Short Put mit Strike=Cap
#
# @params S underlying price
# @params S0 underlying start price
# @params X strike price (bonus level)
# @params B barrier
# @params Cap cap
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params ratio
CappedReverseBonusCertificate<-function(S, S0, X, B, Cap, Time, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  # 1. Short-Position in Aktie : entspricht einem Long-Put mit Strike 2*S0
  ##price1= pmax(2*S0*exp(-T*(r-r_d)) - S,0)
  put <- GBSOption(TypeFlag="p", S, X=2*S0, Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(put,"price"),0)
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. Up-And-Out-Call with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options)
  if (barrierHit)
  {
    price2 <- 0  
  }
  else
  {
    up_out_call <- StandardBarrierOption(TypeFlag="cuo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price2 <- pmax(attr(up_out_call,"price"),0)
    price2 <- ifelse(is.na(price2), 0, price2)
  }
  
  # 3. European-Short-Put with Strike X (=Cap)
  if (Time == 0)
  {
    price3 <- pmax(-1 * (S-Cap), 0)
  }
  else
  {  
    short_put <- GBSOption(TypeFlag="p", S, X=Cap, Time, r, b=r-r_d, sigma) 
    price3 <- pmax(attr(short_put,"price"),0)
  }
  (price1 + price2 - price3) * ratio
}

# Pricing of SprintCertificate (Double Chance)
#
# Duplication:
# (1) Zero-Strike-Call
# (2) 1*Long-Call mit Strike  = X mit Partizipationsrate p
# (3) 2*Short-Call mit Strike = Cap mit Partizipationsrate p
#
# @params S   underlying price
# @params X   strike price
# @params Cap cap
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params participation participation e.g. 1.1 = 110 %
# @params ratio
SprintCertificate<-function(S, X, Cap, Time, r, r_d, sigma, participation, ratio=1)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  
  # 2. Long-Call with Strike X and ratio = (participation - 1)
  long_call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
  price2 <- attr(long_call,"price") * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 3. 2 * Short-Call with Strike Cap and ratio = (participation - 1)
  short_call <- GBSOption(TypeFlag="c", S, X=Cap, Time, r, b=r-r_d, sigma)
  price3 <- attr(short_call,"price") * (participation-1)
  price3 <- ifelse(is.na(price3), 0, price3)
    
  (price1 + price2 - 2 * price3)*ratio
}

# Pricing of OutperformanceCertificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Call
#
# @params S underlying price
# @params X strike price
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params participation participation e.g. 1.1 = 110 %
# @params ratio
OutperformanceCertificate<-function(S, X, Time, r, r_d, sigma, participation, ratio=1)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. Call with Strike X and ratio = (participation - 1)
  call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
  price2 <- attr(call,"price") * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  (price1 + price2)*ratio
}

# Pricing of Capped OutperformanceCertificate
#
# Duplication:
# (1) 1 Zero-Strike-Call
# (2) 1 Call mit Strike X
# (3) 2 Short Call mit Strike Cap
#
# @params S underlying price
# @params X strike price
# @params cap cap
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params participation participation e.g. 1.1 = 110 %
# @params ratio
CappedOutperformanceCertificate<-function(S, X, Cap, Time, r, r_d, sigma, participation, ratio=1)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. Call with Strike X and ratio = (participation - 1)
  call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
  price2 <- attr(call,"price") * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 3. European Short-Call with Strike X (=Cap)
  short_call <- GBSOption(TypeFlag="c", S, X=Cap, Time, r, b=r-r_d, sigma) 
  price3 <- attr(short_call,"price") * (participation-1)
  price3 <- ifelse(is.na(price3), 0, price3)
  
  (price1 + price2 - 2*price3) * ratio
}

# Pricing of OutperformancePlusCertificate
#
# Duplication:
# (1) Zero-Strike-Call
# (2) Call
# (3) Down-And-Out-Put
#
# @params S underlying price
# @params X strike price
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params participation participation e.g. 1.1 = 110 %
# @params ratio
OutperformancePlusCertificate<-function(S, X, B, Time, r, r_d, sigma, participation, ratio=1, barrierHit=FALSE)
{
  # 1. Zero-Strike Call
  zero_strike_call = GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. Call with Strike X and ratio = (participation - 1)
  call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
  price2 <- attr(call,"price") * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 3. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options)
  if (barrierHit)
  {
	  price3 <- 0  
  }
  else
  {
    down_out_put <- StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price3 <- pmax(attr(down_out_put,"price"),0)
    price3 <- ifelse(is.na(price3), 0, price3)
  }
  (price1 + price2 + price3)*ratio
}

# Pricing of Twin-Win Certificate
#
# Duplication:
# (1) Long Zero-Strike-Call in Aktie
# (2) Long Call mit Strike X und BV=(1 - Participation Rate)
# (3) 2xDown-And-Out-Put mit Strike X und Barrier B
#
# @params S underlying price
# @params X strike price (bonus level)
# @params B barrier
# @params Cap cap
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a.
# @params participation
# @params ratio
TwinWinCertificate<-function(S, X, B, Time, r, r_d, sigma, participation = 1, ratio=1)
{
  # 1. Zero-Strike Call
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)  
  
  # 2. Call with Strike X and ratio = (participation - 1)
  call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
  price2 <- attr(call,"price") * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 2. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  down_out_put     = StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
  price3 <- pmax(attr(down_out_put,"price"),0)
  price3 <- ifelse(is.na(price3), 0, price3)
  
  (price1 + price2 + 2 * price3) * ratio
}


# Pricing of Easy-Express-Certificate
#
# Duplication:
# (1) Rückzahlungsbetrag S0 (z.B. 120 EUR)
# (2) Cash-or-Nothing (short-put) mit Strike = B und Auszahlung von (S0-B) : Ist das Underlying am Bewertungstag unter der Barriere, verliert man erstmal noch den Bonus.
# (3) Short-Put mit Strike = B
#
#
# @params S underlying price
# @params S0 Rückzahlungsbetrag bzw. Höchstbetrag (z.B. 120 EUR)
# @params B barrier (B <= S0)
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
EasyExpressCertificate<-function(S, S0, B, Time, r, r_d, sigma, ratio=1)
{
  # 1. Short Cash-or-Nothing-Put mit Strike X=B und Auszahlung (S0-B)
  con_short_put<-CashOrNothingOption("p", S, X=B, (S0-B), Time, r, b=r-r_d, sigma)
  price1=pmax(attr(con_short_put,"price"),0)
  
  # 2. Plain-Vanilla-Short-Put mit Strike B
  plain_short_put<- GBSOption(TypeFlag="p", S, X=B, Time, r, b=r-r_d, sigma)
  price2=pmax(attr(plain_short_put,"price"),0)
  
  (S0 * exp(-r*Time) - price1 - price2) * ratio
}

# Alternative Duplikation von Andy:
#
# (1) Zero-Strike-Call Long
# (2) Cash-or-Nothing Call Long mit Strike B und Auszahlung (S0-B)
# (3) Short-Call mit Strike = B
EasyExpressCertificate2<-function(S, S0, B, Time, r, r_d, sigma, ratio=1)
{
  # 1. Zero-Strike-Call Long
  zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
  price1 <- attr(zero_strike_call,"price")
  price1 <- ifelse(is.na(price1), 0, price1)  
  
  # 2. Cash-Or-Nothing Call Long
  con_long_call <- CashOrNothingOption("c", S, X=B, (S0-B), Time, r, b=r-r_d, sigma)
  price2 <- pmax(attr(con_long_call,"price"),0)
  price2 <- ifelse(is.na(price2), 0, price2)  
  
  # 3. Short-Call mit Strike = B
  plain_short_call<- GBSOption(TypeFlag="c", S, X=B, Time, r, b=r-r_d, sigma)
  price3 <- pmax(attr(plain_short_call,"price"),0)
  price3 <- ifelse(is.na(price3), 0, price3)    
  
  (price1 + price2 - price3) * ratio
}




# Berechnung eines Reverse Convertibles
#
# @param S 
# @param Cap 
# @param Time
# @param r
# @param r_d   
# @param sigma   
# @param nominal 
# @param coupon
# 
ReverseConvertible<-function(S, Cap, Time, r, r_d, sigma, nominal, coupon)
{
  # Anleihenkomponente diskontiert
  price1 <- nominal*(1+coupon*Time)/(1+r)^Time
  
  # Bezugsverhältnis
  ratio <- nominal/Cap
  
  # Short Put
  short_put<- GBSOption(TypeFlag="p", S, X=Cap, Time, r, b=r-r_d, sigma)
  price2 <- attr(short_put,"price")
  price2 <- ifelse(is.na(price2), 0, price2)
  
  (price1 - price2*ratio)/nominal*100
}

# Berechnung eines Reverse Convertibles Plus Pro (RC+PRO)
#
# @param S 
# @param Cap 
# @param B Barrier 
# @param Time
# @param r
# @param r_d   
# @param sigma   
# @param nominal 
# @param coupon
# 
ReverseConvertiblePlusPro<-function(S, Cap, B, Time, r, r_d, sigma, nominal, coupon, barrierHit=FALSE)
{
  # 1. Anleihenkomponente diskontiert
  price1 <- nominal*(1+coupon*Time)/(1+r)^Time
  
  # Bezugsverhältnis
  ratio <- nominal/Cap
  
  # 2. Short Put
  short_put <- GBSOption(TypeFlag="p", S, X=Cap, Time, r, b=r-r_d, sigma)
  price2 <- pmax(attr(short_put,"price"),0)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  # 3. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  if (!barrierHit)
  {
    down_out_put     = StandardBarrierOption(TypeFlag="pdo", S, X=Cap, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price3 <- pmax(attr(down_out_put,"price"),0)
    price3 <- ifelse(is.na(price3), 0, price3)
  }
  else
  {
    price3 <- 0
  }
  
  (price1 - price2*ratio + price3*ratio)/nominal*100
}


# Pricing of Reverse-Outperformance-Certificate
#
# Duplication:
# (1) Short Position in Aktie
# (2) European Short-Put
#
# @param S   underlying price
# @param S0 Startkurs
# @param X   strike price of the embedded option (cap)
# @param Time time to maturity in years
# @param r interest rate p.a. as 0.02 = 2%
# @param r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params participation Partizipationsrate
# @param ratio
ReverseOutperformanceCertificate<-function(S, S0, X, Time, r, r_d, sigma, participation=1, ratio=1)
{
  # 1. Short-Position in Aktie
  put <- GBSOption(TypeFlag="p", S, X=2*S0, Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(put,"price"),0)
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. European Long-Put with Strike X und ratio = (participation-1)
  long_put <- GBSOption(TypeFlag="p", S, X, Time, r, b=r-r_d, sigma) 
  price2 <- pmax(attr(long_put,"price"),0) * (participation-1)
  price2 <- ifelse(is.na(price2), 0, price2)
  
  (price1 + price2) * ratio
}

# Berechnung eines "Airbag-Zertifikats" : Keine Pfadabhängigkeit.
#
# @param S underlying price
# @param X Strike price
# @param B = Absicherungslevel 
# @param Time
# @param r
# @param r_d   
# @param sigma   
# @param participation 
# @param ratio
# 
AirbagCertificate<-function(S, X, B, Time, r, r_d, sigma, participation, ratio=1)
{
  # 1. Cash-Position auf Garantielevel X
  price1 <- X * exp(-r * Time) 
  
  if (Time==0)
  {
    price2 <- pmax(S-X,0)
    price3 <- pmax(B-S,0) * X/B
  }
  else
  {
    # 2. Long Call mit Partizipationsrate für Kurse über Partizipationslevel
    long_call     <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
    price2=attr(long_call,"price") * participation
    
    # 3. Short Put mit Strike B und BV X/B
    short_put     <- GBSOption(TypeFlag="p", S, B, Time, r, b=r-r_d, sigma)
    price3 <- attr(short_put,"price") * X/B
  }
  
  (price1 + price2 - price3) * ratio
}

# Berechnung eines "AirbagPlus-Zertifikats" : Keine Pfadabhängigkeit.
# Societe Generale:
# Ein Airbag Plus-Zertifikat 
# stellt eine Kombination eines Bonus- und eines Airbag-Zertifikates dar. 
# Der Kapitalschutz bleibt selbst dann aktiv, 
# wenn der Index während der Laufzeit die Bonus-Barriere unterschreitet, bei Laufzeitende aber über dieser notiert. Berührt oder unterschreitet der Basiswert während der Laufzeit nie die Barriere, so generiert das Zertifikat die gleichen Auszahlungen wie ein Bonus-Zertifikat. 
# Ein zu Emission festgelegter Cap begrenzt den möglichen Ertrag des Zertifikats.
#
#
# Duplikation:
# 1. Cash-Position auf Garantielevel
# 2. Long Call mit Partizipationsrate für Kurse über Partizipationslevel
# 3. Short Put mit Strike B und BV X/B
# 4. Down-And-Out-Put mit Strike B und unterer Barriere B2
#
# @param S underlying price
# @param X Strike price
# @param B = Absicherungslevel 
# @param Time
# @param r
# @param r_d   
# @param sigma   
# @param participation 
# @param ratio
# @param barrierHit
AirbagPlusCertificate <- function(S, X, B, B2, Time, r, r_d, sigma, participation, ratio=1, barrierHit=FALSE)
{
  # 1. Cash-Position auf Garantielevel
  price1 <- X * exp(-r * Time)
  
  if (Time==0)
  {
	# 2. Long Call mit Partizipationsrate für Kurse über Partizipationslevel
	price2 <- pmax(S-X,0)
	# 3. Short Put mit Strike B und BV X/B
	price3 <- pmax(B-S,0) * X/B
	
	# 4. Down-And-Out-Put with Strike X und Barrier B2
	price4  <- ifelse(B2 <= S & S <= X,  X-S,  0) * X/B
  }
  else
  {
    # 2. Long Call mit Partizipationsrate für Kurse über Partizipationslevel
    long_call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
    price2 <- attr(long_call,"price") * participation
    
    # 3. Short Put mit Strike B und BV X/B
    short_put <- GBSOption(TypeFlag="p", S, B, Time, r, b=r-r_d, sigma)
    price3 <- attr(short_put,"price") * X/B
    
    # 4. Down-And-Out-Put with Strike X und Barrier B2
	  if (!barrierHit)
	  {	
      put <- StandardBarrierOption(TypeFlag="pdo", S, X=X, H=B2, K=0, Time, r, b=r-r_d, sigma)
      price4 <- pmax(attr(put,"price"),0)
    }
    else
	{
	  price4 <- 0
	}
  }
  (price1 + price2 - price3 + price4) * ratio
}

# Berechnung eines "Garantie-Zertifikats" : Garantierte Rückzahlung von einem Betrag.
# 
# Die Anleihe wird am Faelligkeitstag mindestens zu 100% des Nennbetrages zurueckgezahlt (Kapitalgarantie). 
# Uebersteigt der durchschnittliche Schlusskurs der Aktie am Ausuebungstag den Basispreis, 
# so partizipiert der Anleger zusaetzlich on Hoehe der Partizipation an einer positiven Entwicklung der zugrundeliegenden Aktie 
# und erhaelt einen Zusatzbetrag. 
# Dabei wird der durchschnittliche Schlusskurs auf der Basis von der Deutschen Boerse AG im Xetra-Handelssystem offiziell festgestellte 
# Schlusskurse der Aktie an den Bewertungstagen 
# ( 04.04.2006, 04.07.2006, 04.10.2006, 04.01.2007, 04.04.2007, 04.07.2007, 04.10.2007, 04.01.2008, 04.04.2008, 04.07.2008, 06.10.2008, 05.01.2009) berechnet, 
# Der Zusatzbetrag wird wie folgt berechnet : Zusatzbetrag = Partizipation * durchschnittlicher Schlusskurs - Basispreis / Basispreis * Nennbetrag.
#
# @param S underlying price
# @param X strike price
# @param B = Absicherungslevel 
# @param Time
# @param r
# @param r_d   
# @param sigma   
# @param nominal 
# @param coupon
# 
GarantieCertificate<-function(S, X, Time, r, r_d, sigma, participation, ratio=1, nominal)
{
  # 1. Cash-Position auf Garantielevel
  price1 <- nominal * exp(-r * Time)
  
  # 2. Long Call mit Partizipationsrate für Kurse über Basispreis
  if (Time == 0)
  {
    price2 <- pmax(S-X,0) * participation
  }
  else
  {
    long_call <- GBSOption(TypeFlag="c", S, X, Time, r, b=r-r_d, sigma)
    price2 <- attr(long_call,"price") * participation
  }
    
  (price1 + price2) * ratio
}


# Zertifikat zahlt mehrere Boni (Bewertungstage), wenn das Underlying bis zum Bewertungstag nicht unter die Schwelle gerutscht ist.
# --> Summe von Binary-Barrier-Optionen.
# Am Ende der Laufzeit verhält sich das Zertifikat wie ein Discountzertifikat mit Cap.
#
# Returnzertifikate (z.B. DE000CB6VKC2) zahlen einen festen Bonus (in EUR), solange der Basiswert während der Laufzeit die Kursschwelle nicht berührt oder unterschreitet. 
# Nach oben ist das Zertifikat gecappt. Wesentlicher Unterschied zu einem Capped-Bonuszertifikat ist daher, dass der Anleger auch von Kursrückgängen betroffen ist.
# Bei Capped-Bonuszertifikaten wird ein festes Bonusniveau K bezahlt, egal wo die Aktie S steht (sofern B<=S<=X). Damit ist der Aufschlag auf die Aktie unterschiedlich hoch, 
# je nachdem, wo die Aktie innerhalb des Bonusbereiches steht. 
# Hier bei diesem Zertifikat wird ein fester Aufschlag von 11 EUR bezahlt (wenn S>=B während der Laufzeit), plus dem Wert der Aktie von 50, 60 oder 70 EUR.
# Returnzertifikate sind damit nur dann besser als das Direktinvestment, wenn der Basiswert die Kursschwelle nicht berührt und nicht mehr als die Bonuszahlung steigt. 
# Einen absoluten Gewinn erziehlt der Anleger nur, wenn das Underlying nicht fällt. Returnzertifikate zielen damit ganz klar auf Seitwärtsmärkte ab.
#
# Duplikation
#
# (1) Zero-Strike-Call
# (2) Short Call mit Cap
# (3) BinaryBarrierOption mit Auszahlung Bonus, Down-and-Out-Cash-or-Nothing
# 
ReturnCertificate<-function(S, Bonus, S0, B, Cap, Time, r, r_d, sigma, ratio = 1, barrierHit=FALSE)
{
  if (Time == 0)
  {
    if (barrierHit == TRUE)
    {
      price <- pmin(S,Cap)
    }
    else
    {
      price <- pmin(S,Cap) + ifelse(S>B, Bonus[length(Bonus)], 0)
    }
  }
  else
  {
    # 1. Zero-Strike Call
    zero_strike_call <- GBSOption(TypeFlag="c", S, X=0, Time, r, b=r-r_d, sigma) 
    price1 <- attr(zero_strike_call,"price")
    
    # 2. European Short-Call with Strike X (=Cap)
    short_call <- GBSOption(TypeFlag="c", S, X=Cap, Time, r, b=r-r_d, sigma) 
    price2 <- attr(short_call,"price")
    
    # 3. Bonuszahlungen als Binary-Barrier-Optionen
    price3 <- 0
    for (i in 1:length(Time))
    {
      # 3. BinaryBarrierOption mit Auszahlung Bonus : down-and-out-cash-or-nothing
      if (length(Bonus)>1)
      {
        payment <- Bonus[i]
      }
      else
      {
        payment <- Bonus
      }
      binary_cash_or_nothing <- BinaryBarrierOption(TypeFlag = "9", S = S, X = B-0.0001, H = B, K = payment, Time = Time[i], r = r, b = r-r_d, sigma = sigma)
      price3 <- price3 + attr(binary_cash_or_nothing,"price")
    }
    price <- (price1 - price2 + price3)
  }
  price*ratio
}

# Beispiel: RLZ(c("16.06.2008","16.06.2009","16.06.2010"))

# Leveraged Bonus vgl. DESG03DY=SGED
#
LeveragedBonusCertificate<-function(S, X, B, B2, Time, r, r_d, sigma, ratio=1, barrierHit=FALSE)
{
  # 1. Hebelprodukt Call (cdo)
  down_out_call <- StandardBarrierOption(TypeFlag="cdo", S, X=B2, H=B2, K=0, Time, r, b=r-r_d, sigma) 
  price1 <- pmax(attr(down_out_call,"price"),0)
  price1 <- ifelse(is.na(price1), 0, price1)
  
  # 2. Down-And-Out-Put with Strike X and Barrier B (Cash Rebate K = 0 for standard barrier options) 
  if (barrierHit)
  {
    price2 = 0
  }
  else
  {
    down_out_put  <- StandardBarrierOption(TypeFlag="pdo", S, X, H=B, K=0, Time, r, b=r-r_d, sigma) 
    price2 <- pmax(attr(down_out_put,"price"),0)
    price2 <- ifelse(is.na(price2), 0, price2)    
  }
  
  (price1 + price2) * ratio
}

###############################################################################################################
#
# Hilfsfunktionen
#
###############################################################################################################

# Totalverlustwahrscheinlichkeit für Plain-Vanilla-Optionen : Risiko, am Ende der Laufzeit aus dem Geld zu sein. Formal : P(S_T <= X)
#
# @param type : "c" = call, "p" = put
# @param S underlying price
# @param X strike price
# @param T time to maturity 
# @param r riskfree rate
# @param r_d dividend yield
# @param sigma volatility
shortfall_risk<-function(type="c", S, X, T, r, r_d, sigma)
{ 
	mu=r-r_d
	p = pnorm((log(X/S) - (mu - sigma^2/2)*T) / sqrt(sigma^2 * T), mean=0, sd=1);
	if (type == "c") # P(S_T <= X)
	{
		p
	}
	else # P(S_T >= X)
	{
		(1-p)
	}
}

# Hilfsfunktion zum Berechnen der Restlaufzeit in Jahren (RLZ) für eine Reihe von Daten
# Beispiel: RLZ(c("16.06.2008","16.06.2009","16.06.2010"))
#
# @param Vektor mit Datumsangabe
# @param Datumsformat, Standard ist "%d.%m.%Y"
.RLZ<-function(dates, dateformat="%d.%m.%Y", start=NA)
{
	# Referenzdatum t0 ist aktuelles Datum
	t0 = ifelse(!is.na(start) & start != "", as.Date(start, format=dateformat), Sys.Date())
	
	# Restlaufzeiten in Jahren
	rlz = as.numeric(difftime(as.Date(dates ,format=dateformat), as.Date("1970-01-01") + t0), units="days")/365
	rlz
}

annualizeYield<-function(yield, T)
{
	ifelse (T<=1, (yield/T), ((1+yield)^(1/T)-1))
}