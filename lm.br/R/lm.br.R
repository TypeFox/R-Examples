


lm.br  <- function( formula, type ="LL", data, subset,
  weights, inverse =FALSE, var.known =FALSE, na.action,
  contrasts =NULL, offset, ... )  {

## pre-process the input

  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match( c( "formula", "data", "subset", "weights",
    "na.action", "offset" ), names(mf), 0L )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  offset <- as.vector(model.offset(mf))

  type <- toupper(type)

  x <- model.matrix( mt, mf, contrasts )

  y <- model.response( mf, "numeric")
  if( NCOL(y) > 1 )  
    stop( "mutiple response vectors not supported" )
  ynm <- names(mf)[1]

  w <- model.weights(mf)
  w_ <- w
  if( is.vector(w) )  { 
    if(inverse) { 
      if(any(w==0)) 
        stop( "zero variances not allowed" )  
      else
        w <- 1/w 
    }
  }  else
  if( is.matrix(w) )  {
    n <- NROW(y)
    if( any( dim(w) != c(n, n) ) )  
      stop( "dim('weights') invalid" )
    eW <- eigen(w, TRUE)
    d <- eW$values
    if( any(d <= 0) ) stop("'weights' not positive-definite")
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))
    if(!is.null(offset)) 
      offset <- as.vector( A %*% as.matrix(offset) )
  }

# use 'lm.fit' or 'lm.wfit' to check input, check x-rank per 'tolerance',
# and return some of the output list
  z  <-  if( is.null(w) )
      lm.fit( x, y, offset=offset, ... )
    else  {
      if( is.vector(w) )
        lm.wfit( x, y, w, offset=offset,  ... )
      else
        lm.fit( A %*% x, A %*% y, offset=offset, ... )
    }


  if( length(z$coef) > 0  &&  !is.na(z$coef[1]) )  {

# 'x1' is the term with a coefficient changepoint
    dn <- colnames(x)
    xint <- if( dn[1] == "(Intercept)" )  TRUE  else  FALSE
    if(xint) x1c=2 else x1c=1
    x1 <- as.vector( x[,x1c] )
    x1nm <- dn[x1c]


# 'xb' is the 'x' input matrix, but with two vectors for 'x1'
#  before and after changepoint,  and with zeroes for columns
#  that were found to be linearly dependent
    nx <- ncol(x)
    nxb <- ncol(x)+1
    xb <- matrix( 0, nrow(x), nxb )
    xb[,1] <- x[,1]
    for(i in x1c:nx) 
      xb[,i+1] <- if( is.na(z$coef[i]) )  0  else  x[,i]

    bnm <- paste( " ", x1nm, "< theta" )
    bpnm <- paste( " ", x1nm, "> theta" )
    rownames(xb) <- rownames(x)
    colnames(xb) <-  if( xint )
      { if( x1c < nx )  
          c("alpha", bnm, bpnm, dn[(x1c+1):nx])  
        else  
          c("alpha", bnm, bpnm) 
      }
      else
      { if( x1c < nx )  
          c(bnm, bpnm, dn[(x1c+1):nx])  
        else
          c(bnm, bpnm) 
      }


#  re-order for 'x1' non-decreasing,  drop rows with  w = 0,  
#  drop columns that are linearly dependent
    x_ <- x[ order(x1), , drop = FALSE]
    y_ <- y
    if(!is.null(offset)) y_ <- y_ - as.vector(model.offset(mf))
    y_ <- y_[ order(x1) ]
    if( !is.null(w_) )  
      w_ <- if(is.matrix(w_))  w_[ order(x1), order(x1) ]  
            else  w_[ order(x1) ]

    if( is.vector(w_) )  if( any(w_==0) )  {
      ok <- w_!=0
      w_ <- w_[ok]
      y_ <- y_[ok]
      x_ <- x_[ok, , drop = FALSE]
    }
    if( is.vector(w_) )  w_ <- as.matrix(w_)
    if( is.null(w_) )  w_ <- as.matrix( -1 )

    x_fullcol <- x_
    x_ <- x_[ , !is.na(z$coef), drop = FALSE ]

    if( type=="LL" ) {
      if( xint )  model_num <- 1  else
        stop( "'alpha'=0 not supported for type \"LL\"" )
    } else {
      if( type=="TL" ) {
        if( xint )  model_num <- 2  else  model_num <- 3
      }  else  {
        if( type=="LT" ) {
          if( xint )  model_num <- -2  else  model_num <- -3 
        }  else  
          stop( "'type' must be \"LL\", \"LT\" or \"TL\"" )
      }
    }


#  Drop columns that are dependent at some changepoint value, via 'while' loop.
#  If the x-matrix is dependent at  changepoint = 'th'  then  
#  Q*f(th)=0  on lower rows,  where  f(th) = max( x1-th, 0 ) and QR = x-matrix.
#  The C++ object includes a subroutine that returns 'th' for minimum Q*f(th).
    x_dep <- TRUE

    while( x_dep )  {

#  construct the C++ object
      obj <- new( Cpp_Clmbr, y_, x_, w_, model_num, 
            as.integer(inverse), as.integer(var.known) )

      thQfmin <- obj$param()[5]
      xb[ ,x1c] <- if( type=='TL' )  0  else
          { if( is.infinite(thQfmin) )  -1  else  pmin(x1 - thQfmin, 0 ) }
      xb[ ,x1c+1] <- if( type=='LT' )  0  else
          { if ( is.infinite(thQfmin) )  1  else  pmax(x1 - thQfmin, 0 ) }

      z  <-  if( is.null(w) )
          lm.fit( xb, y, offset=offset, ... )
        else  {
          if( is.vector(w) )
            lm.wfit( xb, y, w, offset=offset,  ... )
          else
            lm.fit( A %*% xb, A %*% y, offset=offset, ... )
        }

      xbr <- z$rank
      nb <- if(type=='LL')  ncol(x_) + 1  else  ncol(x_)
      if( xbr < nb )  {
        xb[ , is.na(z$coef) ]  <- 0
        ok <- vector( 'logical', nx )
        ok[] <- TRUE
        for(i in (x1c+1):nx)  if( is.na(z$coef[i+1]) )  ok[i]<- FALSE
        x_ <- x_fullcol[ , ok, drop = FALSE ]
      } else
        x_dep <- FALSE

    }


# output
    par <- obj$param()
    if( is.na( par[1] ) )  {   # case of perfect line
      xb <- x
      if(xint) colnames(xb)[1] <- "  alpha"
      for(i in 1:nx) if( is.na(z$coef[i+1]) ) xb[,i] <- 0
    }  
    else  {
      xb[ , x1c ]  <-  if( type=='TL' )  0  else
          { if( is.infinite(par[1]) )  1  else  pmin( x1-par[1], 0 ) }
      xb[ , x1c+1 ]  <-  if( type=='LT' )  0  else
          { if ( is.infinite(par[1]) )  1  else  pmax( x1-par[1], 0 ) }
    }

    if( is.infinite(par[1]) )  {  # could occur for models with alpha=0
      if( type=='LT' )
        { colnames(xb)[1] <- "  (Intercept)";  colnames(xb)[2] <- bnm }
      if( type=='TL' )
        { colnames(xb)[1] <- bpnm;  colnames(xb)[2] <- "  (Intercept)" }
    }


    z  <-  if( is.null(w) )
        lm.fit( xb, y, offset=offset,  ... )
      else  {
        if( is.vector(w) )
          lm.wfit( xb, y, w, offset=offset,  ... )
        else
          lm.fit( A %*% xb, A %*% y, offset=offset, ... )
      }


    if( type=='TL' )  z$coefficients[x1c] <- 0
    if( type=='LT' )  z$coefficients[x1c+1] <- 0

#  add 'theta' to the coefficients
    for(i in ncol(xb):1 )  {
      z$coef[i+1] <- z$coef[i]
      names(z$coef)[i+1] <- names(z$coef)[i]
    }
    z$coef[1] <- par[1]
    names(z$coef)[1] <- "theta"


#  accessor functions
    z$CppObj  <-  obj

    z$ci <- function(...)  .ci( z, ... )

    z$cr <- function(...)  .cr( z, ... )

    z$sl <- function(...)  .sl( z, ... ) 

    z$mle <- function( )  (z$CppObj)$mle( )

    z$sety <- function( rWy )  .sety( z, rWy )


#  output list
    p <-  if( is.na( par[1] ) )  z$rank  else  z$rank + 1
    z$no_of_parameters <- p
    z$df.residual <- nrow(x_) - p
    z$xint <- xint
    z$x1 <- x1
    z$x1nm <- x1nm
  }

  if( is.vector(w) )  {
    if(!is.null(offset)) z$rWoffset <- sqrt(w) * offset
  }  else
  if( is.matrix(w) )  {
    z$fitted.values <- drop(Ainv %*% z$fitted.values)
    z$residuals <- drop(Ainv %*% z$residuals)
    if(!is.null(offset)) z$rWoffset <- offset
  }

  class(z) <- "lm.br"
  z$call <- call
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$terms <- mt
  z$xlevels <- .getXlevels(mt, mf)
  z$offset <- as.vector(model.offset(mf))
  z$y <- y
  z$ynm <- ynm
  z$weights <- model.weights(mf)
  z$type <- type
  z
}




