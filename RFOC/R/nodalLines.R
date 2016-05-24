nodalLines <-
function( strike,dip,rake, PLOT=TRUE)
{
                   
### *   function plots nodal lines in the
### *      lower-hemisphere equal-area projection           
### *                                         
### ************************************
###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
  ### translated to R by Jonathan M. Lees
  if(missing(PLOT)) { PLOT=TRUE }
  
 N = length(strike)
  
n.1 =matrix(ncol=3, nrow=N)
u.1 =matrix(ncol=3, nrow=N)
n.1[,1] = -sin(dip*pi/180)*sin(strike*pi/180)
n.1[,2] =  sin(dip*pi/180)*cos(strike*pi/180)
n.1[,3] = -cos(dip*pi/180)

u.1[,1] =  cos(rake*pi/180)*cos(strike*pi/180)+
  cos(dip*pi/180)*sin(rake*pi/180)*sin(strike*pi/180)
 
u.1[,2] =  cos(rake*pi/180)*sin(strike*pi/180)-
  cos(dip*pi/180)*sin(rake*pi/180)*cos(strike*pi/180)
 
u.1[,3] = -sin(rake*pi/180)*sin(dip*pi/180)

N = length(strike)

### ----------------------------------------------
###  lower hemisphere equal-area projection
### -------------------------------------
projekce = -1  
### ------------------------------------
###  1st nodal lines
### -----------------------------------
ksi.min = 0.1
ksi.max = 360.1
ksi.step = 1
xyang = seq(from=ksi.min, by=ksi.step, to=ksi.max)
nang = length(xyang)
x = rep(NA, nang )
y= rep(NA, nang )

sqrt2 = sqrt(2.)
n = n.1

for(j  in  1:N)
  {
    n1 = sqrt(n[j,1]^2+n[j,2]^2);
    n3 = n[j,3];
    m1 = n[j,1]/n1;
    m3 = n[j,2]/n1
    if (n3<0){ n3 = -n3; m1 = -m1; m3 = -m3    }
###  the vertical component must be always negative!
    
    k = 1
    for (  ksi in xyang )
      {
        k1 = sin(ksi*pi/180);
        k3 = cos(ksi*pi/180)

        smer = -c( - k1*m1*n3 + k3*m3, -(k1*m3*n3 + k3*m1),k1*n1 )
        
        if(smer[3]<0)
          { ###  plot of one hemispehre only
            theta = acos(k1*n1)*180/pi
            if (theta<1.e-6)
              { ###  theta must be non-zero
                fi = abs(atan(smer[1]/smer[,2])*180/pi)
                sin.fi = sign(smer[1])*sin(fi)
                cos.fi = sign(smer[2])*cos(fi)
              }
            else
              {
                sin.fi = smer[1]/sin(theta*pi/180)
                cos.fi = smer[2]/sin(theta*pi/180)
              } #            
            x[k] = sqrt2*projekce*sin(theta*pi/360)*cos.fi
            y[k] = sqrt2*projekce*sin(theta*pi/360)*sin.fi
            
           # if (k>1)
           #   {
           #     lines(c(x[k-1], x[k]),c(y[k-1], y[k]) ) 
           #   } #
            k=k+1
          }
        else
          {
            k=1 
          }
      }
    if(k>1)
      {
        PLANE1 = list(x=x, y=y)

      }
  }

### -----------------------------------------
###  2nd nodal lines
### -----------------------------------------
n = u.1

for(j  in  1:N)
  {
    n1 = sqrt(n[j,1]^2+n[j,2]^2);
    n3 = n[j,3];
    m1 = n[j,1]/n1 ;
    m3 = n[j,2]/n1 
    if(n3<0) { n3 = -n3; m1 = -m1; m3 = -m3    }
###  the vertical component must be always negative!
    
    k = 1
    for (  ksi in xyang )
      {
        k1 = sin(ksi*pi/180)
        k3 = cos(ksi*pi/180)

        smer = -c(- k1*m1*n3 + k3*m3, -(k1*m3*n3 + k3*m1),k1*n1)
        if(smer[3]<0)
          { ###  plot of one hemisphere only
            theta = acos(k1*n1)*180/pi
            if (theta<1.e-6)
              {###  theta must be non-zero
                fi = abs(atan(smer[1]/smer[,2])*180/pi)
                sin.fi = sign(smer[1])*sin(fi)
                cos.fi = sign(smer[2])*cos(fi)
              }
            else
              {
                sin.fi = smer[1]/sin(theta*pi/180)
                cos.fi = smer[2]/sin(theta*pi/180)
               }             
            x[k] =sqrt2*projekce*sin(theta*pi/360)*cos.fi
            y[k] =sqrt2*projekce*sin(theta*pi/360)*sin.fi
            
          #  if (k>1)
           #   {
           #     lines(c(x[k-1], x[k]),c(y[k-1], y[k]) ) 
          #     }
            k=k+1
          }
        else
          {
            k=1 
           }
       }

    if(k>1)
      {
        PLANE2 = list(x=x, y=y)
        
      }
   }


return(list(PLANE1=PLANE1, PLANE2=PLANE2) )

}
