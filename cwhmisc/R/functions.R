inrange <- function(x,r) (!is.na(x) && min(r) <= x & x <= max(r))

mod <- function( m, n  ) m %% n

modS <- function( m, n ) return  (m - (trunc( m/n )*n))

modR <- function(x, y ) return  (x - (floor( x/y )*y))

zero <- function( x ) 0.0

one  <- function( x ) 1.0

onebyx <- function( x ) 1.0/x
  
sqr <- function( x ) x^2

powr <- function( a, x ) a^x

equal <- function( x, y )  x == y

equalFuzzy <- function( x, y, prec=8*.Machine$double.eps, rel=TRUE ) abs(x - y ) <= prec*(ifelse(rel,abs(x)+abs(y),1.0))

quotmean <- function(x,y) mean(x,na.rm = TRUE)/mean(y,na.rm = TRUE)

safeDiv <- function( num, den ) {
  q <- ifelse (num==0 & den==0, 1,  num/den)
  return( ifelse (is.infinite(q), c3Q, q) )
}

solveQeq <- function(a,b,c) { # solve ax^2 + bx + c = 0 for x
  ab <- abs(b) + abs(c)
  if (Mod(a)*0.125 + ab == ab) { 
    # a almost 0 to within 8*.Machine$double.eps
    res <- -c/b
    if (is.infinite(res)) res <- c(NA,NA) 
  } else {
    dis <- b^2 - 4*a*c
    if ((Im(dis) == 0) && (Re(dis) < 0.0)) {  # complex case
      dis <- complex(real=dis)
    }
    x1 <- -signp(b)*(abs(b) + sqrt(dis))/a*0.5
    res <- if (dis!=0.0) c(x1,c/a/x1) else c(x1,x1)
  }
  return (res)
}

chsvd <- function(s){ s$u %*% diag(s$d) %*% t(s$v) }

divmod <- function(x,y) c(x %/% y, x %% y)

submod <- function(x, v) {
  if (x <= 0) return ( c(0,0) )
  v <- c(0, v)
  ii <- sum( v < x  )
  return( c( ii-1, x - v[ii]) )
}  # submod

dsm <-  function( x, w ) {
   dm <- divmod( x, sum(w) )
   sm <- submod( dm[2], cumsum(w) )
   res <- c( dm[1], sm )
   return( res )
 } # dsm

