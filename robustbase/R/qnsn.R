### Note: till 2010, a slightly wrong constant = 2.2219 was in use.
### Error reported by Peter Ruckdeschel, U.Bayreuth, 15.Juli 2010
### correct constant == 1 / (sqrt(2) * qnorm(5/8)) == 2.219144
### -- but wrong constant, 2.2219, is already in the the original Fortran qn.f
Qn.corr <- 2.2219 / 2.21914
##' Qn finite sample correction factor  (not exported, but "available")
##' Version 1
Qn.finite.c <- function(n)
    (if (n %% 2) 1.6069 +(-2.333 - 3.1/n)/n # n  odd
     else        3.6667 +( 2.591 - 5.5/n)/n # n  even
    )/n + 1
## Version built on res <- cbind(Res.sml, Res.mid) ## and the models there
Qn.finite.c <- function(n)
    (if (n %% 2) 1.60188 +(-2.1284 - 5.172/n)/n # n  odd
     else        3.67561 +( 1.9654 +(6.987 - 77/n)/n)/n # n  even
    )/n + 1

Qn <- function(x, constant = 2.21914, finite.corr = missing(constant))
{
    ## Purpose: Rousseeuw and Croux's  Q_n() robust scale estimator
    ## Author: Martin Maechler, Date: 14 Mar 2002, 10:43
    n <- length(x)
    if(n == 0) return(NA) else if(n == 1) return(0.)

    r <- constant *
        .C(Qn0, as.double(x), n, res = double(1))$res

    if (finite.corr) {
	if (n <= 12) ## n in 2:12 --> n-1 in 1:11
	    ## n=2: E[Q_2] = E|X - Y| = sqrt(pi)/2, fc = sqrt(pi)/2/2.21914
	    r* c(.399356, # ~= fc = 0.3993560233
                 ## These are from MM's simulation("Res3"), Nsim = 2^27 ~= 134 mio:
                 ## ~/R/MM/Pkg-ex/robustbase/Qn-simulation.R
                 .99365, .51321, .84401, .61220,
                 .85877, .66993, .87344, .72014,
                 .88906, .75743)[n - 1L]
	else
	    r / Qn.finite.c(n)
    }
    else r
}

## This is the old version -- available for back "compatibility":
Qn.old <- function(x, constant = 2.2219, finite.corr = missing(constant))
{
    ## Purpose: Rousseeuw and Croux's  Q_n() robust scale estimator
    ## Author: Martin Maechler, Date: 14 Mar 2002, 10:43
    n <- length(x)
    if(n == 0) return(NA) else if(n == 1) return(0.)

    r <- constant *
        .C(Qn0, as.double(x), n, res = double(1))$res

    if (finite.corr)
	(if (n <= 9) { # n in 2:9 --> n-1 in 1:8
            c(.399,.994, .512,.844, .611,.857, .669,.872)[n - 1]
	} else {
	    if (n %% 2) ## n odd
		n / (n + 1.4)
	    else ## n even
		n / (n + 3.8)
        }
         ) * r
    else r
}


Sn <- function(x, constant = 1.1926, finite.corr = missing(constant))
{
    ## Purpose: Rousseeuw and Croux's  S_n() robust scale estimator
    ## Author: Martin Maechler, Date: 14 Mar 2002, 10:43

    n <- length(x)
    if(n == 0) return(NA) else if(n == 1) return(0.)

    r <- constant * .C(Sn0,
                       as.double(x), n,
                       as.integer(!is.unsorted(x)),# is.sorted
                       res = double(1), a2 = double(n))$res
    ## NB: a2[] could be used for confidence intervals and other estimates!
    if (finite.corr) (
	if (n <= 9) {
            c(0.743, # n = 2
              1.851, 0.954,# 3 & 4
              1.351, 0.993,# 5 & 6
              1.198, 1.005,# 7 & 8
              1.131)[n - 1]
	} else if (n %% 2) # n odd, >= 11
            n / (n - 0.9)
        else # n even, >= 10
            1
    ) * r
    else r
}

wgt.himedian <- function(x, weights = rep(1,n))
{
    ## Purpose: weighted hiMedian of x
    ## Author: Martin Maechler, Date: 14 Mar 2002
    n <- length(x <- as.double(x))
    stopifnot(storage.mode(weights) %in% c("integer", "double"))
    if(n != length(weights))
	stop("'weights' must have same length as 'x'")
    ## if(is.integer(weights)) message("using integer weights")
    .C(if(is.integer(weights)) wgt_himed_i else wgt_himed,
       x, n, weights,
       res = double(1))$res
}

## To be used directly as  'scaleFun'  in  'covOGK()' :
s_Qn <- function(x, mu.too = FALSE, ...)
    c(if(mu.too) median(x), Qn(x, ...))

s_Sn <- function(x, mu.too = FALSE, ...)
    c(if(mu.too) median(x), Sn(x, ...))
