Hd <- function( h, m, s ) {
  return ( s/60.0 + m )/60.0 + h
} ## end  Hd

Hms <- function( hd ) {
  h <- int( hd )
  s <- hd - h;  s <- s*60.0
  m <- int( s )
  s <- 60.0 * (s - m)
  return (list(h=h,m=m,s=s))      
} ## end  Hms

Hdms <- function( hd ) {
  r <- Hms( hd )
  return ((r$s/100.0 + r$m)/100.0 + r$h)
} ## end  Hdms

Hmsd <- function( hms ) {
  h <- int( hms )
  s <- hms - h;  s <- s*100.0
  m <- int( s )
  s <- 100.0 * (s - m)
  return ((s/60.0 + m )/60.0 + h)
} ## end  Hmsd
