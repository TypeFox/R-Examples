`zebra` <-
function(x,y,Dx, dx, dy, lab="", pos=1, col=c("black", "white"), cex=1, textcol="black", xpd=TRUE,  PMAT=NULL)
{
###  zebra: draw a zebra bar km scale on a map; Km ; KM scale; kmscale
###  x, y = location of left point 
###   DX  =length km
###    dx  = alternating distance
###   dy = thickness of zebra


  if(missing(PMAT)) { PMAT=NULL  }
  if(missing(lab)) { lab=NULL  }
  if(missing(col)) { col=c("black", "white")  }
  if(missing(cex)) { cex=1  }
  if(missing(pos)) { pos=1 }
  if(missing(xpd)) { xpd=TRUE }
  if(missing(textcol)) { textcol="black"  }

  
  for(i in 1:(Dx/dx))
    {
      px1 = x+(i-1)*dx
      py1 = y
      px2 = x+(i)*dx
      py2 = y+dy

      if(!is.null(PMAT))
        {
          tem1 = trans3d(px1, py1, rep(0, length(py1)) , PMAT)
          tem2 = trans3d(px2, py2, rep(0, length(py2)) , PMAT)
          rect(tem1$x[1], tem1$y[1], tem2$x[1], tem2$y[1], col=col[i%%2 + 1] , xpd=xpd)

        }
      else
        {
          rect(px1, py1, px2, py2, col=col[i%%2 + 1], xpd=xpd)
        }
      
    }

  ytext = y
  if(pos==3)
    {
      ytext = y+dy
      
    }
      
  if(!is.null(PMAT))
    {

      tem1 = trans3d(x, ytext, rep(0, length(y)) , PMAT)
      text(tem1$x[1], tem1$y[1], labels="0", pos=pos, cex=cex, xpd=xpd, col=textcol)
      
      tem1 = trans3d(px2, ytext, rep(0, length(y)) , PMAT)
      text(tem1$x[1], tem1$y[1], labels=Dx, pos=pos, cex=cex, xpd=xpd, col=textcol )
      
      tem1 = trans3d((x+px2)/2, ytext, rep(0, length(y)) , PMAT)
      text(tem1$x[1], tem1$y[1], labels=lab, pos=pos, cex=cex , xpd=xpd, col=textcol)
      
    }
  else
    {
      text(x,ytext, labels="0", pos=pos, cex=cex, xpd=xpd, col=textcol)
      text(px2, ytext, labels=Dx, pos=pos, cex=cex, xpd=xpd, col=textcol )
      text((x+px2)/2, ytext, labels=lab, pos=pos, cex=cex, xpd=xpd, col=textcol)
    }

}

