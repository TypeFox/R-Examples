IntInNode <- function(layout,cex,shape,m,width=0.2,triangles=TRUE,col="black",side=1,inside=TRUE)
{
  N <- nrow(layout)
  if (length(cex)==1) cex <- rep(cex,N)
  if (length(shape)==1) shape <- rep(shape,N)
  if (length(col)==1) col <- rep(col,N)
  if (length(side)==1) side <- rep(side,N)
  
  # m is vector of margins to plot lines, NA indicates no line
  # side: 1. bottom, 2. left, 3. top, 4. right.
  # inside: if TRUE thresholds are plotted in the node, filling from top to bottom, if FALSE they are plotted at the side.
  
  for (i in seq_along(m))
  {
    if (!is.na(m[i]))
    {
      #       browser()
      x <- layout[i,1]
      y <- layout[i,2]
      xran <- qgraph:::Cent2Edge(layout[i,1],layout[i,2],pi/2,cex[i],cex[i],shape[i])[1] - x
      yran <- qgraph:::Cent2Edge(layout[i,1],layout[i,2],0,cex[i],cex[i],shape[i])[2] - y
      
      if (!inside)
      {
        if (side[i]==1)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran+m[[i]][j]*xran*2,x-xran+m[[i]][j]*xran*2),c(y-yran-width*yran,y-yran+width*yran),col=col[i])
          }
        } else if (side[i]==2)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran-width*xran,x-xran+width*xran),c(y-yran+m[[i]][j]*yran*2,y-yran+m[[i]][j]*yran*2),col=col[i])
          }        
        } else if (side[i]==3)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran+m[[i]][j]*xran*2,x-xran+m[[i]][j]*xran*2),c(y+yran-width*yran,y+yran+width*yran),col=col[i])
          }        
        } else if (side[i]==4)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x+xran-width*xran,x+xran+width*xran),c(y-yran+m[[i]][j]*yran*2,y-yran+m[[i]][j]*yran*2),col=col[i])
          }        
        }
      } else 
      {
        if (side[i]==1)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran+m[[i]][j]*xran*2,x-xran+m[[i]][j]*xran*2),c(y-yran,y+yran),col=col[i])
          }
        } else if (side[i]==2)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran,x+xran),c(y-yran+m[[i]][j]*yran*2,y-yran+m[[i]][j]*yran*2),col=col[i])
          }        
        } else if (side[i]==3)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran+m[[i]][j]*xran*2,x-xran+m[[i]][j]*xran*2),c(y-yran,y+yran),col=col[i])
          }        
        } else if (side[i]==4)
        {
          for (j in 1:length(m[[i]]))
          {
            lines(c(x-xran,x+xran),c(y-yran+m[[i]][j]*yran*2,y-yran+m[[i]][j]*yran*2),col=col[i])
          }        
        }        
      }
    }
  }
}
#           if (triangles)
#           {
#             points(x,y-yran+m[[i]][j]*yran*2,pch=17,cex=cex[1]/10,col=col[i])
#           }