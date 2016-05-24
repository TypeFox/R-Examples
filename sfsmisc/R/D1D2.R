## This is also sym.linked into
## Martin's WpDensity package /u/maechler/R/Pkgs/WpDensity/

###------- Numerical Derivatives ------------------------------------------

### Test Programs and examples for those two are in
### -->  "/u/maechler/S/NUMERICS/D1-tst.S"
###
### For 'optimal' 2nd Deriv.:  d2.est(..)
###  --> "/u/maechler/S/NUMERICS/diff2.S"  "/u/maechler/S/NUMERICS/diff2-user.S"


D1tr <- function(y, x = 1)
{
  ## Purpose:  discrete trivial estimate of 1st derivative.
  ## -------------------------------------------------------------------------
  ## Arguments: x is optional
  ## -------------------------------------------------------------------------
  ##--> See also D1.naive in ~/S/D1-tst.S (and the (smoothing) one: 'D1') !
  ## Author: Martin Maechler, ~ 1990
  n <- length(y)
  if(length(x) == 1)
    c(y[2] - y[1], 0.5 * (y[-(1:2)] - y[-((n-1):n)]), y[n] - y[n-1])/x
  else {
    if(n != length(x)) stop("lengths of 'x' & 'y' must equal")
    if(is.unsorted(x)) stop("'x' must be sorted !")
    c(y[2] - y[1], 0.5 * (y[-(1:2)] - y[-((n-1):n)]), y[n] - y[n-1]) /
      c(x[2] - x[1], 0.5 * (x[-(1:2)] - x[-((n-1):n)]), x[n] - x[n-1])
  }
}


D1ss <- function(x, y, xout = x, spar.offset = 0.1384, spl.spar=NULL)
{
  ## Purpose: Numerical first derivatives of  f() for   y_i = f(x_i) + e_i.
  ## Find  f'(xout)  -- using smoothing splines with GCV'
  ## Author: Martin Maechler, Date:  6 Sep 92, 00:04
  ## -------------------------------------------------------------------------
  ## Arguments: x = { x_i } MUST be sorted increasingly // y = { y_i }
  ## -------------------------------------------------------------------------
  sp <-
    if(is.null(spl.spar)) {
      sp <- smooth.spline(x,y)
      smooth.spline(x,y, spar = sp$ spar + spar.offset)
    } else smooth.spline(x,y, spar = spl.spar)
  predict(sp, xout, deriv = 1) $ y
}

D2ss <- function(x, y, xout = x, spar.offset = 0.1384, spl.spar=NULL)
{
  ## Purpose: Numerical 2nd derivative of  f() for   y_i = f(x_i) + e_i.
  ##          Find  f''(xout) -- using smoothing splines (with GCV) -- DOUBLY:
  ##          f --ss-> f' --ss-> f''
  ## -------------------------------------------------------------------------
  ## Arguments: x = { x_i } MUST be sorted increasingly // y = { y_i }
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 29 Jan 1997, 17:55 -- for S-plus
  ## -------------------------------------------------------------------------

  use.fudge <- is.null(spl.spar)
  if(use.fudge) { ##-- use  GCV * 'spar.offset' ---
    if(is.null(spar.offset)) stop("must specify 'spl.spar' OR 'spar.offset'!")
    lf <- length(spar.offset)
    if(!is.numeric(spar.offset) || lf == 0 || lf > 2)
      stop("'spar.offset' must be numeric(1 or 2) !")
    if(lf == 1) spar.offset <- rep(spar.offset, 2)
    sp <- smooth.spline(x,y)
    sp <- smooth.spline(x,y, spar = spar.offset[1] + sp $ spar)
    spl.spar <- numeric(2); spl.spar[1] <- sp $ spar
  }
  else {
    lf <- length(spl.spar)
    if(!is.numeric(spl.spar) || lf == 0 || lf > 2)
      stop("'spl.spar' must be numeric(1 or 2) !")
    if(lf == 1) spl.spar <- rep(spl.spar, 2)
    sp <- smooth.spline(x,y, spar = spl.spar[1])
  }

  D1 <- predict(sp, x, deriv = 1) $ y #-- 1st derivative ...

  if(use.fudge) { ##-- use  GCV * 'spar.offset' ---
    sp <- smooth.spline(x, D1)
    sp <- smooth.spline(x, D1, spar = spar.offset[2] + sp $ spar)
    spl.spar[2] <- sp $ spar
  } else {
    sp <- smooth.spline(x, D1, spar = spl.spar[2])
  }
  if(is.unsorted(xout))
        xout <- sort(xout)
  list(x=xout, y = predict(sp, xout, deriv = 1) $ y,
       spl.spar = spl.spar, spar.offset = spar.offset)
}



D1D2 <- function(x, y, xout = x, spar.offset = 0.1384,
                 deriv = 1:2, spl.spar=NULL)
{
    ## Purpose: Numerical first derivatives of  f() for   y_i = f(x_i) + e_i.
    ## Find  f'(xout) & f''(xout) -- using smoothing splines with GCV'
    ## Author: Martin Maechler, Date:  23 Sep 1992, 9:40ith GCV'
    ## Author: Martin Maechler, Date:  23 Sep 1992, 9:40
    ## -------------------------------------------------------------------------
    ## Arguments: x = { x_i } MUST be sorted increasingly // y = { y_i }
    ## -------------------------------------------------------------------------

    if(is.unsorted(xout))
        xout <- sort(xout)
    sp <-
        if(is.null(spl.spar)) {
            sp <- smooth.spline(x,y)
            smooth.spline(x,y, spar = sp$ spar + spar.offset)
        } else smooth.spline(x,y, spar = spl.spar)
    c(list(x = xout,
           D1 = if(any(deriv==1)) predict(sp, xout, deriv = 1) $ y,
           D2 = if(any(deriv==2)) predict(sp, xout, deriv = 2) $ y),
      sp[c("spar", "df")])
}
