getreachind <- function(net,x,y, ...) UseMethod("getreachind")


# distance of point (x,y) from reach (x1,y1) to (x2,y2)

get.dist <- function(x,y,x1,y1,x2,y2)
{
  if ( is.na(x) | is.na(y) | is.na(x1) | is.na(y1) | is.na(x2) | is.na(y2) ) return(NA)
  
  # distance between points:
  
  if ( x1==x2 & y1==y2 ) return(sqrt((x-x1)^2+(y-y1)^2))
  
  # distance to reach:
  
  # calculate t with minimal distance for xline = x1 + t*(x2-x1), yline = y1 + t*(y2-y1)
  
  t <- ( (x-x1)*(x2-x1) + (y-y1)*(y2-y1) ) / ( (x2-x1)^2 + (y2-y1)^2 )
  
  if ( t < 0 ) return(sqrt((x-x1)^2+(y-y1)^2))  # distance to (x1,y1)
  if ( t > 1 ) return(sqrt((x-x2)^2+(y-y2)^2))  # distance to (x2,y2)
  return(sqrt((x1+t*(x2-x1)-x)^2+(y1+t*(y2-y1)-y)^2))
}


getreachind.rivernet <- function(net,x,y,...)
{
  n <- length(x)
  if ( n != length(y) ) return(NA)
  if ( n == 0 )         return(NA)
  
  ind  <- rep(NA,n)
  dist <- rep(NA,n)
  for ( i in 1:n )
  {
    for ( j in 1:length(net$reaches) )
    {
      for ( k in 1:(length(net$reaches[[j]]$x)-1) )
      {
        dist.new <- get.dist(x[i],y[i],
                             net$reaches[[j]]$x[k]  ,net$reaches[[j]]$y[k],
                             net$reaches[[j]]$x[k+1],net$reaches[[j]]$y[k+1])
        if ( !is.na(dist.new) )
        {
          if ( is.na(dist[i]) )
          {
            dist[i] <- dist.new
            ind[i]  <- j
          }  
          else
          {
            if ( dist.new < dist[i] )
            {
              dist[i] <- dist.new
              ind[i]  <- j
            }
          }
        }
      }
    }
  }
  return(data.frame(reach.ind=ind,reach.dist=dist))
}
