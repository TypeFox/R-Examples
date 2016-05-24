PTaxes <-
function( strike,dip,rake)
{

### *   function plots P and T axes in the lower-hemisphere equal-area projection
### *                                                 
### *************************************************
###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
### translated to R by Jonathan M. Lees
  
  Vnorm <-function(X){ return( sqrt(sum(X^2)) ) }

N = length(strike)
  
n.1 =matrix(ncol=3, nrow=N)
u.1 =matrix(ncol=3, nrow=N)
P.osa=matrix(ncol=3, nrow=N)
T.osa=matrix(ncol=3, nrow=N)
P.azimuth= rep(NA, N)
T.azimuth= rep(NA, N)
T.theta= rep(NA, N)
P.theta= rep(NA, N)
P.x= rep(NA, N)
P.y= rep(NA, N)
T.x= rep(NA, N)
T.y= rep(NA, N)

n.1[,1] = -sin(dip*pi/180)*sin(strike*pi/180)
n.1[,2] =  sin(dip*pi/180)*cos(strike*pi/180)
n.1[,3] = -cos(dip*pi/180)

u.1[,1] =  cos(rake*pi/180)*cos(strike*pi/180) +
  cos(dip*pi/180)*sin(rake*pi/180)*sin(strike*pi/180)

u.1[,2] =  cos(rake*pi/180)*sin(strike*pi/180) -
  cos(dip*pi/180)*sin(rake*pi/180)*cos(strike*pi/180)

u.1[,3] = -sin(rake*pi/180)*sin(dip*pi/180)

### --------------------------------------------
###  lower hemisphere equal-area projection
### --------------------------------------------
projekce = -1  

### ---------------------------------------------
###  P/T axes
### -------------------------------------------

for (  i in 1 : N )
  {
    P.osa[i,] = (n.1[i,]-u.1[i,])/Vnorm(n.1[i,]-u.1[i,])
    T.osa[i,] = (n.1[i,]+u.1[i,])/Vnorm(n.1[i,]+u.1[i,])
    
    if (P.osa[i,3]>0) { P.osa[i,1]=-P.osa[i,1];
                        P.osa[i,2]=-P.osa[i,2];
                        P.osa[i,3]=-P.osa[i,3]    }
    if (T.osa[i,3]>0) {  T.osa[i,1]=-T.osa[i,1];
                         T.osa[i,2]=-T.osa[i,2];
                         T.osa[i,3]=-T.osa[i,3]    }

    fi = atan(abs(P.osa[i,1]/P.osa[i,2]))*180/pi

    if (P.osa[i,1]>0 & P.osa[i,2]>0) {  P.azimuth[i] = fi      } #  1. kvadrant
    if (P.osa[i,1]>0 & P.osa[i,2]<0)  { P.azimuth[i] = 180-fi  } #  2. kvadrant
    if (P.osa[i,1]<0 & P.osa[i,2]<0)  { P.azimuth[i] = fi+180  } #  3. kvadrant
    if (P.osa[i,1]<0 & P.osa[i,2]>0)  {  P.azimuth[i] = 360-fi } #  4. kvadrant

    P.theta[i] = acos(abs(P.osa[i,3]))*180/pi

    fi = atan(abs(T.osa[i,1]/T.osa[i,2]))*180/pi

    if (T.osa[i,1]>0 & T.osa[i,2]>0){  T.azimuth[i] = fi   } #  1. kvadrant
    if (T.osa[i,1]>0 & T.osa[i,2]<0){ T.azimuth[i] = 180-fi } #  2. kvadrant
    if (T.osa[i,1]<0 & T.osa[i,2]<0){  T.azimuth[i] = fi+180 } #  3. kvadrant
    if (T.osa[i,1]<0 & T.osa[i,2]>0){  T.azimuth[i] = 360-fi } #  4. kvadrant

    T.theta[i] = acos(abs(T.osa[i,3]))*180/pi

    P.x[i] = sqrt(2.)*projekce*sin(P.theta[i]*pi/360)*sin(P.azimuth[i]*pi/180)
    P.y[i] = sqrt(2.)*projekce*sin(P.theta[i]*pi/360)*cos(P.azimuth[i]*pi/180)

    T.x[i] = sqrt(2.)*projekce*sin(T.theta[i]*pi/360)*sin(T.azimuth[i]*pi/180)
    T.y[i] = sqrt(2.)*projekce*sin(T.theta[i]*pi/360)*cos(T.azimuth[i]*pi/180)

   } #

points(P.y,P.x,col='green')
text(P.y,P.x,labels='P', pos=3)
points(T.y,T.x,col='green', pch=3)
text(T.y,T.x,labels='T', pos=4)


return(NULL)

}
