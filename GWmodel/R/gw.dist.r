################################################################################
# gw.dist: caculate a distance vector or distance matrix between (a) regression point and data points
# dp.locat: numeric vector of two colunms, coordinates of data points
# rp.locat: numeric vector of two colunms, coordinates of regression points
# focus: an integer, indexing to the current regression point, if focus=0, all the distances between all the regression points and data points will be calculated and a distance matrix will be returned;
# if 0<focus<length(rp.locat), then the distances between the focusth regression points and data points will be calculated and a distance vector will be returned;
# p: the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
# theta: an angle in radian to roate the coordinate system, default is 0
# longlat: if TRUE, great circle will be caculated

 gw.dist<- function(dp.locat, rp.locat, focus=0, p=2, theta=0, longlat=F)
 {
   if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
      stop("Please input correct coordinates of data points")
  
   if (!missing(rp.locat)) 
     {
       rp.given<-T
     }
    else
     {
       rp.given<-F 
       rp.locat<- dp.locat
     } 
   
   if (!is.numeric(rp.locat)||dim(rp.locat)[2]!=2)
      stop("Please input correct coordinates of regression points")
   if (focus<0||focus>length(rp.locat[,1]))
      stop("No regression point is fixed")

   n.rp<-length(rp.locat[,1])
   n.dp<-length(dp.locat[,1])
   if (focus>0)
       dist.res<-numeric(n.dp)
   else
       dist.res<-matrix(numeric(n.rp*n.dp),nrow=n.dp)
   ###Rotate the coordiante system
   if (p!=2||theta!=0||!longlat)
   {
     dp.locat<-coordinate.rotation(dp.locat, theta)
     rp.locat<-coordinate.rotation(rp.locat, theta)
   }
   ####Calculate a distance vector at regression point "focus"
   if (focus>0)
   {
      if (longlat)   ###Great circle distance, and FORTRAN or C++ code to be added
         dist.res<-spDistsN1(dp.locat, rp.locat[focus,],longlat=longlat)
      else
      { 
          for (i in 1:n.dp)
          {
            
            if (p==2)
               dist.res[i]<-dist(rbind(dp.locat[i,],rp.locat[focus,]))
            else if (is.infinite(p))
               dist.res[i]<-dist(rbind(dp.locat[i,],rp.locat[focus,]),method="maximum")
            else
               dist.res[i]<-dist(rbind(dp.locat[i,],rp.locat[focus,]),method="minkowski",p=p)
          }
      }               
   }
   else
   {
      if(rp.given) ####Regression points are given, which means different sets of points are used for fitting and sampling
      {
        for (i in 1:n.rp)
        {
          if (longlat)
             dist.res[,i]<-spDistsN1(dp.locat, matrix(rp.locat[i,],nrow=1),longlat=longlat)
          else
          {
              for (j in 1:n.dp)
              {
                if (p==2)
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]))
                else if (is.infinite(p))
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]),method="maximum")
                else
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]),method="minkowski",p=p)
              }
          }
        }
      }
      else
      {
        for (i in 1:(n.rp-1))
        {
          for (j in (i+1):n.dp)
          {
            if (longlat)
            {
               dist.res[j,i]<-spDistsN1(matrix(dp.locat[j,],nrow=1), matrix(rp.locat[i,],nrow=1),longlat=longlat)
            }
            else
            {
                if  (p==2)
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]))
                else if (is.infinite(p))
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]),method="maximum")
                else
                   dist.res[j,i]<-dist(rbind(dp.locat[j,],rp.locat[i,]),method="minkowski",p=p)
            }
            dist.res[i,j]<-dist.res[j,i]
          }
        }
      }
   }
        
   
   dist.res     

 }
 
 
 
 ##Transform coordinates after rotating the coordinate axes by theta
 ##This function should be writted in FORTRAN or C++
coordinate.rotation<-function(coords, theta)
{
  n<-as.integer(nrow(coords))
  Xlist<-as.double(coords[,1])
  Ylist<-as.double(coords[,2])
  rotated.x<-Xlist*cos(theta)-Ylist*sin(theta)
  rotated.y<-Xlist*sin(theta)+Ylist*cos(theta)
  rotated.xy<-cbind(rotated.x,rotated.y)
  rotated.xy
}

####Chebyshev distance
Chebyshev<-function(dp.locat, rp.focus)
{
  N<-length(dp.locat[,1])
  dist.V<-numeric(N)
  for (i in 1:N)
   dist.V[i]<-max(c(abs(dp.locat[i,1]-rp.focus[1]),abs(dp.locat[i,2]-rp.focus[2])))

  dist.V
}

