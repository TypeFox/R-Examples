blf2 <-
function(A,b,c,delta,l,u)
{
 ##  require(bvls)
  ZEE = bvls::bvls(A,b,l,u)
  x0=ZEE$x
  
  gamma0=t(c) %*% x0
  alpha=1.0e-3
  
###  First, make sure that || A x0 - b || <= delta.  If not, then the problem
###  is just infeasible.
  
  vtest = Vnorm(A %*% x0 - b)
  
  if(vtest > delta)
    {
      print(paste("problem is just infeasible", vtest , " > ", delta ))
      
      xmin=NULL
      xmax=NULL
      return(list( xmin=xmin,xmax=xmax ) )
    } 
  
  
###  now find xmax
###  Find a gamma so that || A*x(gamma)-b || > delta
  step=0.1*abs(gamma0)+0.01
  gamma=gamma0+step
  G=rbind( A ,  alpha %*% t(c) )
  
  d=c(b , alpha*gamma)
   ZEE = bvls::bvls(G,d,l,u)
  xgamma=ZEE$x
  
  xgammanorm=Vnorm( A %*% xgamma - b)

###  left is a step added to gamma known to have xgammanorm < delta
  left=0
  while((xgammanorm < delta) & (step < 1.0e10))
    {
      left = step
      step=step*2
      gamma=gamma0+step
      G=rbind( A ,  alpha %*% t(c) )
      
      d=c(b,  alpha%*%gamma)

      ZEE = bvls::bvls(G,d,l,u)
      xgamma=ZEE$x
      
      xgammanorm=Vnorm(A %*% xgamma - b)
    }




  right=step

###   Now, do a binary search to narrow the range.
  while((right-left)/(0.0001+right) >0.01)
    {
      mid=(left+right)/2
      gamma=gamma0+mid
      G=rbind(A,  alpha%*% t(c) ) 
      d=c( b,  alpha %*% gamma)

      ZEE = bvls::bvls(G,d,l,u)
      xgamma=ZEE$x
      xgammanorm=Vnorm(A %*% xgamma-b)
      if(xgammanorm > delta)
        {
          right=mid
        }
      else
        {
          left=mid
        } 
    }
  
###  Save this solution as xmax
  gamma=gamma0+left
  G=rbind(A , alpha %*% t(c)) 
  d=c(b,  alpha%*% gamma)
  ZEE  =   bvls::bvls(G,d,l,u)
  xmax=ZEE$x
  cmax=t(c) %*% xmax


###  Now find xmin
###  Find a gamma so that || A*x(gamma)-b || > delta
### 
  step=0.1*abs(gamma0)+0.01
  gamma=gamma0-step
  G=rbind(A,  alpha%*% t(c) )
  d=c( b,  alpha %*% gamma )
  ZEE = bvls::bvls(G,d,l,u)
  xgamma=ZEE$x
  xgammanorm=Vnorm(A %*% xgamma-b)

###  left is a step added to gamma known to have xgammanorm < delta
  left=0
  while((xgammanorm < delta) & (step < 1.0e10))
    {
      left = step
      step=step*2
      gamma=gamma0-step
      G=rbind(A,  alpha %*% t(c) )
      d=  c(b,  alpha %*% gamma)

      ZEE = bvls::bvls(G,d,l,u)
      xgamma=ZEE$x
      xgammanorm=Vnorm(A %*% xgamma-b)
    } 
  right=step

###   Now, do a binary search to narrow the range.
  while((right-left)/(0.0001+right) >0.01)
    {
      mid=(left+right)/2
      gamma=gamma0-mid
      G= rbind(A ,alpha %*% t(c) )
      d=c(b,  alpha %*% gamma)

      ZEE =bvls::bvls(G,d,l,u)
      xgamma= ZEE$x
      xgammanorm=Vnorm(A %*% xgamma - b )
      if (xgammanorm > delta)
        {
          right=mid
        }
      else
        {
          left=mid
        } 
    } 

###  Save this solution as xmin
  gamma=gamma0-left
  G=rbind( A,  alpha %*% t(c) )
  d=c(b,  alpha %*% gamma)

  ZEE = bvls::bvls(G,d,l,u)
  xmin=ZEE$x
  cmin=t(c) %*% xmin

  
###  ensure that xmin corresponded with the lower dot average
  if(cmin > cmax)
    {
      temp=xmax
      xmax=xmin
      xmin=temp
    } #


  return(list( xmin=xmin,xmax=xmax ))


}
