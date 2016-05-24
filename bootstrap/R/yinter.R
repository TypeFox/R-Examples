"yinter"<-
function(x, y, z, increasing = TRUE)
{
    ## for function defined by (x(i),y(i)), i=1,...n, interpolate y
    ## value at z
    ##
    ## used by boott
    
    if(increasing == FALSE) {
        x <- -1 * x
        x <- x[length(x):1]
        y <- y[length(y):1]
        z <- -1 * z
    }
    n <- length(x)
    if (z <= x[1]) {
        ii <- 1;jj <- 2;while ( y[jj]==y[ii]) {jj <- jj+1;}
    }
    else if (z>=x[n]) {
        jj <- n;ii <- n-1;while ( y[ii]==y[jj]) {ii <- ii-1;}
    }
    else {ii <- 1;
          while  (!((x[ii] <= z) & (z <= x[ii+1])))
          {ii <- ii+1;}
          jj <- ii+1;
      }
    if (x[ii]==x[jj]) {result <- (y[ii]+y[jj])/2;}
    else {slope <- (y[jj]-y[ii])/(x[jj]-x[ii]);
          result <- y[ii]+slope*(z-x[ii]);
      }
    
    return(result)
}
