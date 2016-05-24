makeneighborsw <- function(coords,method="neighbor",m=1,d,cum=TRUE)
{
 # initialisation
 if(ncol(coords)!=2) stop("coords should have 2 columns")
  xc<-coords[,1]
  yc<-coords[,2]

 # condition
  if(length(xc)!=length(yc))
  stop("Number of coords not equal")

  findneighbors<- function (xc,yc,m)
    {
      # last modified 16/09/08
      n <- length(xc)
      nnlist <- matrix(0, nrow=n, ncol=m)
      d <- as.matrix(dist(cbind(xc,yc)))

      # calcul de la matrice contenant les indices des plus proches voisins de chaque observation

      for (i in 1:n)
      {
        d1 <- d[i,-i]
        names(d1) <- c(1:n)[which(c(1:n)!=i)]
        d2 <- sort(d1,index.return=TRUE)
       
        x <- d2$x
        ix <- d2$ix    
    
        ind <- which(x[2:(length(x)-1)]==x[3:length(x)])
        p <- length(ind)
    
        if(p!=0)
        {for (j in 1:p)
          ix[which(x==x[ind[j]+1])]=ix[sample(which(x==x[ind[j]+1]))]
        }
        
        nnlist[i,] <- names(d1[ix[1:m]])
      }
    
    return(nnlist)
    } 
    
  # initialisation
  if(m>length(xc)) m<-length(xc)

  if(method=="distance"||method=="both")
  { W.dist<-dist(cbind(xc,yc))
    W.dist[which(W.dist <= d,arr.ind=TRUE)]<-1
    W.dist[which(W.dist > d,arr.ind=TRUE)]<-0

    W.dist<-as.matrix(W.dist)
  }
 
  if(method=="neighbor"||method=="both")
  {n.list <- findneighbors (xc,yc,m)
   w.neigh <- matrix(0,nrow=length(xc),ncol=length(xc))

     if(cum==TRUE)
     {
      for (i in 1:length(xc))
        {
          w.neigh[i,as.integer(n.list[i,])] <- 1
        }
     }
     else
     {
      for (i in 1:length(xc))
        {
          w.neigh[i,as.integer(n.list[i,m])] <- 1
        }
     } 
  }
 
 
if(method=="neighbor") return(as.matrix(w.neigh)) 
if(method=="distance") return(W.dist)
if(method=="both") return(W.dist * as.matrix(w.neigh))

}

