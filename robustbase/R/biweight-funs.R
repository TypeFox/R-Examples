#### These Chi() and Psi() used to be called by lmrob() functions
####  but no longer --> Have interface  via .psi2ipsi() and .C(..)
#### FIXME: integrate these with the psi-rho objects --> ./psi-rho-funs.R

## In the vignette ../vignettes/psi_functions.Rnw, we call this
##  scaled \rho  "\tilde{\rho}"
##- Maronna et al (2006) define their rho to be standardized
##-    (only if possible <==>  only if redescending psi !)

##- {TODO: *Where* in the Hampel_et_al book ??? }
## Hampel et al (1986):  \chi(x) := \rho(x) / \rho(\infty)
##                       ======
## <==> chi() is a scaled version of rho(.) such that
##  \chi(\infty) = \max_x \chi(x) = 1

## ==> Chi'() is just a scaled version of psi() :
## with current scale (new for psi()):
##	 i)  Chi'(x, c) == (6/c^2) Psi(x,c)
## ==>	 ii) Chi''(x,c) == (6/c^2) Psi'(x,c)
## and       Chi (x, c) == (6/c^2) Rho(x,c), where Psi(.) = Rho'(.)

tukeyChi <- function(x, cc, deriv = 0)
{
    .Deprecated("Mchi")
    x <- x / cc
    x2 <- x*x
    out <- x2 > 1
    switch(deriv + 1,
       {  ## deriv = 0
	   r <- x2*(3 + x2*(-3 + x2))
	   r[out] <- 1
       },
       {  ## deriv = 1
	   r <- 6/cc * x * (1-x2)^2
	   r[out] <- 0
       },
       {  ## deriv = 2
	   r <- 6/(cc^2) * (1 - x2) * (1 - 5*x2)
	   r[out] <- 0
       },
       stop("deriv must be in {0,1,2}"))
    r
}

## we call this  '*Psi1'  such as to not be confounded with
## the (future!) S4 object tukeyPsi() !
tukeyPsi1 <- function(x, cc, deriv = 0)
{
    .Deprecated("Mpsi")
    ## This version of psi() is scaled such that psi'(0) = 1
    x2 <- (x / cc)^2
    if(deriv < 0) out <- x2 > 1 else in. <- x2 < 1
    switch(deriv + 2,
       {  ## deriv = -1
	   c. <- cc^2/6
	   r <- c.*(1 - (1- x2)^3)
	   r[out] <- c.
	   r
       },
       {  ## deriv = 0
	   in. * x * (1-x2)^2
       },
       {  ## deriv = 1
	   in. * (1 - x2) * (1 - 5*x2)
       },
       {  ## deriv = 2
	   in. * 4*x/cc^2 * (5*x2 - 3)
       },
       stop("deriv must be in {-1,0,1,2}"))
}

if(FALSE)
tukeyPsi1Ex <- function (x, cc, deriv = 0)
## tukeyPsi1Ex <- function (x, cc = 4.685, deriv = 0)
##                               ^^^^^^^^^
{
  ## This version of psi() is scaled such that psi'(0) = 1
  u <- pmin((x/cc)^2, 1)
  if(deriv < 0)
    return((1 - (1-u)^3)*cc^2/6)
  if(deriv == 0)
    return(x * (1 - u)^2)
  return((1 - u) * (1 - 5 * u))
}
