`SynAnticline` <-
function(x,y, syn=TRUE, spacing=NULL, N=1, r1= 1, r2= 1.2, h1= 0, h2= 0,endtol=.1,
                    REV=FALSE, col='black', ...)
{
  if(missing(spacing))  spacing=NULL
  if(missing(REV)) REV=FALSE
  if(missing(r1)) { r1 = 1 }
  if(missing(r2)) { r2 = 1.2 }
  if(missing(h1)) { h1 = .5 }
  if(missing(h2)) { h2 = .5 }

  if(missing(col)) { col='black' }
  if(missing(N)) { N = 1 }
    if(missing(syn)) { syn=TRUE }

  
  if(REV){ x= rev(x); y = rev(y) }
     if(missing(endtol)) { endtol=.1 }


  n = length(x)
  
  g = PointsAlong(x, y, N=N,  endtol=endtol)

  

  lines(x,y, col=col, ...)

  arrows(x[n-1], y[n-1], x[n], y[n], col=col, length = 0.1 )
  
 ## g$rot$sn = -g$rot$sn

  if(syn)
    {
      cs  = g$rot$cs
      sn  = -g$rot$sn
    }
  else
    {
      cs  = -g$rot$cs
      sn  = g$rot$sn

    }
  
  g$rot$cs  = sn
  g$rot$sn  = cs
  
  
  horseshoe(g$x  , g$y , r1=r1, r2=r2, h2=h2, h1=h1, rot=g$rot, col=col)
  


  
}

