"xinter"<-
function(x, y, z, increasing = TRUE)
{
    ## for function defined by (x(i),y(i)), i=1,...n, interpolate x
    ## value at z
    if(increasing == FALSE) {
        x <- -1 * x
        x <- x[length(x):1]
        y <- y[length(y):1]
    }
    
    n <- length(x)

    if (z <= y[1]) {
        ii <- 1;jj <- 2;while(x[jj]==x[ii]) {jj <- jj+1;}}
    else if (z >= y[n]) {
        jj <- n;ii <- n-1;while(x[ii]==x[jj]) {ii <- ii-1;}}
    else {
        ii <- 1;                                                             
        while( (!((y[ii] <= z) & (z <= y[ii+1])))                
              &
              (!((y[ii]>= z) & (z>= y[ii+1])))  )   
        {ii <- ii+1;}                                             
        jj <- ii+1;                                                           
    }                                                                   
    if (x[ii]==x[jj]) {result <- (x[ii])}  else 
    if ((y[ii]==y[jj]) & (z <= y[1]))
    {result <- x[1];} else 
    if ((y[ii]==y[jj]) & (z >= y[n]))
    {result <- x[n];} else 
     if (y[ii]==y[jj]) {result <- (x[ii]+x[jj])/2;} else 
         {slope <- (y[jj]-y[ii])/(x[jj]-x[ii]);
          result <- x[ii]+(z-y[ii])/slope;
      }
    
    if(increasing == FALSE) {
        result <- -1 * result
    }
    return(result)
}
