"calcEbspline" <- 
function (theta, data) 
{
  dim(theta)<-c(length(theta)/data@ncomp, data@ncomp)
 
  #knots<-c(seq(660,800,by=10) 


  #knots<-c(545.7008,560.4046,573.1123,585.8131,598.5238,611.2246, 
  #623.9354,636.6362,649.3408,660.0477,670.7554,680.4592)

  #knots<-c(545, 555, 565, 575, 585, 590, 610, 620, 630, 640, 
  #         650, 665)

  bs(data@x2,df=15) %*% theta

}
