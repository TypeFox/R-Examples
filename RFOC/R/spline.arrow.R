spline.arrow<-function(x,y=0,kdiv=20, arrow=1, length=.2, col="black",
                       thick = 0.01, headlength = 0.2, headthick = 0.1, code=2, ...)
{
  ####  make a spline interpolated line from a set of clicks
  ####   add an arrow at the end

  ###   arrow=2 -> fancy arrow
 ### require(GEOmap)
  if(arrow>1)
    {
      ###   should be plotted with asp=1
      ###require(RFOC)
    }
  
  if(is.list(x))
    {
      y = x$y
      x = x$x
    }

  G = GEOmap::getspline(x, y, kdiv=kdiv)
  lines(G, col = col, ...)
  n = length(G$x)
   if(arrow==1)
    {
      arrows(G$x[n-1], G$y[n-1], G$x[n], G$y[n], length=length, col = col)
      if(code==3)
        {
          arrows(G$x[2], G$y[2], G$x[1], G$y[1], length=length, col = col)

        }

      
}
  if(arrow==2)
    {
      
      
      fancyarrows(G$x[n-1], G$y[n-1], G$x[n], G$y[n],
                  thick =thick ,
                  headlength =  headlength,
                  headthick =headthick,
                  col = col , border = col )
    }

  
  invisible(G)
}


