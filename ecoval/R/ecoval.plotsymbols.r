ecoval.plotsymbol <- function(x,y,r,col,square=F)
{
  n.points.seg <- 20
  m <- length(col)
  if ( square & (m==4))
  {
    polygon(x+r*c(0,0,1,1,0)  ,y+r*c(0,1,1,0,0)  ,col=col[1]) 
    polygon(x+r*c(0,0,1,1,0)  ,y+r*c(0,-1,-1,0,0),col=col[2]) 
    polygon(x+r*c(0,0,-1,-1,0),y+r*c(0,-1,-1,0,0),col=col[3]) 
    polygon(x+r*c(0,0,-1,-1,0),y+r*c(0,1,1,0,0)  ,col=col[4]) 
  }
  else
  {
    for ( j in 1:m )
    {
      beta1 <- (j-1)*2*pi/m
      beta2 <-   j  *2*pi/m
      beta  <- beta1 + (0:n.points.seg)/n.points.seg*(beta2-beta1)
      x.seg <- x + c(0,r*sin(beta),0)
      y.seg <- y + c(0,r*cos(beta),0)
      polygon(x.seg,y.seg,col=col[j])
    }
  }
}


ecoval.plotsymbols <- function(nodes,x,y,r,u,
                               square     = F,
                               labels     = NA,
                               col        = utility.calc.colors(),
                               pos.legend = NA,
                               cex.nodes  = 1)
{
  points(x,y,pch=19)
  
  n.points.seg <- 20
  n <- length(x)
  if ( n < 1 ) return
  if ( length(y) != n ) return
  if ( nrow(u) != n ) return
  m <- length(nodes)
  if ( m < 1 ) return
  for ( i in 1:n )
  {
    # determine nearest neighbours:
    
    alpha <- 0.5*pi
    if ( n > 1 )
    {
      ind.others <- (1:n)[-i]
      dist <- rep(NA,n-1)
      for ( j in 1:(n-1) )
      {
        dist[j] <- sqrt( (x[ind.others[j]]-x[i])^2 + (y[ind.others[j]]-y[i])^2 )
      }
      ind.order <- order(dist)
      
      if ( dist[ind.order[1]] < 5*r )
      {
        dx <- x[ind.others[ind.order[1]]] - x[i]
        dy <- y[ind.others[ind.order[1]]] - y[i]
        alpha <- atan2(dx,dy) + pi  # opposite of nearest point
        if ( n > 2 )
        {
          if ( dist[ind.order[2]] < 5*r )
          {
            dx <- x[ind.others[ind.order[2]]] - x[i]
            dy <- y[ind.others[ind.order[2]]] - y[i]
            alpha2 <- atan2(dx,dy)+pi
            if ( abs(alpha-alpha2) > pi ) alpha2 <- alpha2 + 2*pi
            alpha <- 0.5*(alpha+alpha2)  # opposite of mean of nearest points
          }
        }
      }
    }
    delta.x <- 1.5*r*sin(alpha)
    delta.y <- 1.5*r*cos(alpha)
    
    ind <- match(nodes,colnames(u))
    ind.not.na <- !is.na(ind)
    v <- rep(NA,length(nodes))
    v[ind.not.na] <- as.numeric(u[i,ind[ind.not.na]])
    colors <- utility.get.colors(v,col)
    ecoval.plotsymbol(x=x[i]+delta.x,y=y[i]+delta.y,
                   r=r,col=colors,square=square)
    if ( length(labels) >= i )
    {
      if ( !is.na(labels[i]) )
      {
        x.lab <- x[i]+delta.x+r
        pos   <- 4
        if ( delta.x < 0 )
        {
          x.lab <- x[i]+delta.x-r
          pos   <- 2
        }
        y.lab <- y[i]+delta.y
        text(x.lab,y.lab,labels[i],pos=pos,adj=c(0,0.5),
             cex=cex.nodes)
      }
    }
  }
  if ( !is.na(pos.legend[1]) )
  {
    ecoval.plotsymbol(x=pos.legend[1],y=pos.legend[2],r=r,
                   col=rep("white",m),square=square)
    for ( i in 1:m )
    {
      angle <- 2*pi*(i-0.5)/m
      x <- pos.legend[1] + 1.75*r*sin(angle)
      y <- pos.legend[2] + 1.75*r*cos(angle)
      pos <- 4
      if ( x < pos.legend[1] ) pos <- 2
      text(x=x,y=y,labels=nodes[i],pos=pos,
           offset=0,adj=c(0,0.5),cex=cex.nodes)
    }
  }
}

