# $Id: wapply.R 1012 2006-11-14 22:25:06Z ggorjan $

"wapply" <- function( x, y, fun=mean, method="range",
                    width, n=50, drop.na=TRUE, pts, ...)
{
  method <- match.arg(method, c("width","range","nobs","fraction"))
  if(missing(width))
    if( method=="nobs" ) width <- max(5, length(x)/10 )
  else
    width <- 1/10

  if(method == "width" || method == "range" )
    {
      if(method=="range")
        width <- width * diff(range(x))

      if(missing(pts))
        pts <- seq(min(x),max(x),length.out=n)
      
      result <- sapply( pts, function(pts,y,width,fun,XX,...)
                      {
                        low <- min((pts-width/2),max(XX))
                        high <- max((pts+width/2), min(XX))
                        return (fun(y[(XX>= low) & (XX<=high)],...))
                      },
                      y=y,
                      width=width,
                      fun=fun,
                      XX = x,
                      ...)
      if(drop.na)
        {
          missing <- is.na(pts) & is.na(result)
          pts <- pts[!missing]
          result <- result[!missing]
        }
      
      return(list(x=pts,y=result))
    }
  else # method=="nobs" || method=="fraction"
    {
      if( method=="fraction")
        width <- floor(length(x) * width)

      ord <- order(x)
      x  <- x[ord]
      y  <- y[ord]

      n  <- length(x)
      center  <- 1:n
      below <- sapply(center - width/2, function(XX) max(1,XX) )
      above <- sapply(center + width/2, function(XX) min(n,XX) )

      retval  <- list()
      retval$x  <- x
      retval$y  <- apply(cbind(below,above), 1,
                         function(x) fun(y[x[1]:x[2]],...) )
                         
      if(drop.na)
        {
          missing <- is.na(retval$x) | is.na(retval$y)
          retval$x <- retval$x[!missing]
          retval$y <- retval$y[!missing]
        }


      return(retval)
    }
      
}


