# $Id: colorpanel.R 736 2005-11-18 00:16:28Z warnes $

colorpanel <- function(n,low,mid,high)
  {
    if(missing(mid) || missing(high) )
      {
        ## convert to rgb
        low <- col2rgb(low)
        if(missing(high))
          high <- col2rgb(mid)
        else
          high <- col2rgb(high)
        
        red    <- seq(low[1,1], high[1,1], length=n)/255
        green  <- seq(low[3,1], high[3,1], length=n)/255
        blue   <- seq(low[2,1], high[2,1], length=n)/255
      }
    else # use a center color
      {
        isodd <- odd(n)
        if(isodd)
          {
            n <- n+1
          }

        ## convert to rgb
        low <- col2rgb(low)
        mid <- col2rgb(mid)
        high <- col2rgb(high)

        ## determine length of each component
        lower <- floor(n/2)
        upper <- n - lower
        
        red  <- c(
                  seq(low[1,1], mid [1,1], length=lower),
                  seq(mid[1,1], high[1,1], length=upper)
                  )/255
        
        green <- c(
                   seq(low[3,1], mid [3,1], length=lower),
                   seq(mid[3,1], high[3,1], length=upper)
                   )/255
        
        blue <- c(
                  seq(low[2,1], mid [2,1], length=lower),
                  seq(mid[2,1], high[2,1], length=upper)
                  )/255
        
        if(isodd)
          {
            red   <- red  [-(lower+1)]
            green <- green[-(lower+1)]
            blue  <- blue [-(lower+1)]
          }
      }
      
    rgb(red,blue,green)
  }



# Generate red-to-green colorscale
redgreen <- function(n) colorpanel(n, 'red', 'black', 'green')
greenred <- function(n) colorpanel(n, 'green', 'black', 'red' )
bluered  <- function(n) colorpanel(n, 'blue','white','red')
redblue  <- function(n) colorpanel(n, 'red','white','blue')
