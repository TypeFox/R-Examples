######### R-function:dfltCounts  ######### 
# Obtain default set of grid counts from a 
# multivariate point cloud 'x'.
# Last changed: 18 JUL 2005

dfltCounts <- function(x,gridsize=rep(64,NCOL(x)),h=rep(0,NCOL(x)), supp=3.7, range.x, w)
{
   x <- as.matrix(x)
   d <- ncol(x)
   n <- nrow(x)
   if (missing(w)) w <- rep(1,n)
   
   if (missing(range.x))
   {
     range.x <- list()
     for (id in 1:d)
       range.x[[id]] <- c(min(x[,id])-supp*h[id],max(x[,id])+supp*h[id])  
   }

   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))

   gpoints <- list()
   for (id in 1:d)
      gpoints[[id]] <- seq(a[id],b[id],length=gridsize[id])  
 
   if ((d!=1)&(d!=2)&(d!=3)&(d!=4)) stop("binning implemented only for d=1,2,3,4")

   if (d==1) gcounts <- binning(x=x, bgridsize=gridsize, h=h, xmin=a, xmax=b, w=w)$counts
   else if (d>1) gcounts <- binning(x=x, bgridsize=gridsize, H=diag(h^2), xmin=a, xmax=b, w=w)$counts
  
   ##if (d==1) gcounts <- linbin.ks(x,gpoints[[1]], w=w) 
   ##if (d==2) gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]], w=w)
   ##if (d==3) gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]], w=w)
   ##if (d==4) gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]], w=w)
   
   return(list(counts=gcounts,range.x=range.x))
}

######## End of dfltCounts ########
