#Function to evaluate the smoothed fitted density at given points
'dslcd' <- function(x, lcd, A = hatA(lcd)){
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  d <- ncol(lcd$x)
  if(!(is.matrix(A) && (all(eigen(A)$values > .Machine$double.eps^0.5*4)) && (dim(A)[1]==d) && (dim(A)[2]==d))) {
    stop("error: Hat matrix A must be ",d," by ",d," positive definite")
  }
  if(prod(eigen(A)$values) < 1e-7) {
    warning("warning: eigen values of hat matrix A are too small, numerically unstable")
  }
  if(is.vector(x) && length(x)==d) {
    x = matrix(x,ncol=d)
  }
  if( is.vector( x ) && d==1 ) {
    x = matrix( x )
  }

  if(is.matrix(x) && ncol(x)==d)
  {
    if (d==1) val = dslcd_1D (x, lcd, A)
    else if(d==2) val = dslcd_2D (x, lcd, A)
    else val = dslcd_com (x, lcd, A)
    return(val)
  }
  else stop("error: x must be a vector, or a numeric matrix with ",d," columns")
}


'simplex_unit_05_2_nd_pt' <- function( n )
{
# Return the points to evaluate and their weights by Combinatorial Methods
# Modified from the code written by John Burkardt  
# Reference:
#    Axel Grundmann, H M Moller,
#    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
#    SIAM Journal on Numerical Analysis,
#    Volume 15, Number 2, 1978, 282-290.


  x = rep(0.0,n)
# Group 1
  for ( i in 1:n ) x[i] = 1.0 / ( n + 1 )
  coef = ( n + 1 )^4 / ( 32 * ( n + 2 ) * ( n + 3 ))
  ev = c(x,1-sum(x))

# Group 2
  a = 1.0 / ( n + 3 )
  b = 3.0 / ( n + 3 )
  for ( i in 1:n ) x[i] = a
  coef = c(coef,-( n + 3 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 4 ) ))
  ev = rbind(ev, c(x,1-sum(x)))
  for ( i in 1:n )
  {
    x[i] = b
    coef = c(coef,-( n + 3 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 4 ) ))
    ev = rbind(ev, c(x,1-sum(x)))
    x[i] = a
  }

# Group 3
  a = 1.0 / ( n + 5 )
  b = 5.0 / ( n + 5 )
  for ( i in 1:n ) x[i] = a
  coef = c(coef, ( n + 5 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) ))
  ev = rbind(ev, c(x,1-sum(x)))
  for ( i in 1:n )
  {
    x[i] = b
    coef = c(coef, ( n + 5 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) ))
    ev = rbind(ev, c(x,1-sum(x)))
    x[i] = a
  }

# Group 4
  a = 1.0 / ( n + 5 )
  b = 3.0 / ( n + 5 )

  for ( i in 1:n )
  {
    for ( j in 1:n )
    {
      x[j] = a
    }
    x[i] = b
    coef = c(coef, ( n + 5 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) ))
    ev = rbind(ev, c(x,1-sum(x)))
    if ( i+1 <= n)
    {
      for ( j in (i+1):n )
      {
        x[j] = b
        coef = c(coef, ( n + 5 )^4 / ( 16 * ( n + 1 ) * ( n + 2 ) * ( n + 3 ) * ( n + 4 ) ))
        ev = rbind(ev, c(x,1-sum(x)))
        x[j] = a
      }
    }
  }
  ev = matrix(ev,ncol=n+1)
  coef = coef / factorial(n)
  return (list(pt=ev,coef=coef))
}


'dslcd_com' <- function (y, lcd, A) 
{
#  Evaluate the density by Combinatorial Methods
    triang = lcd$triang
    x = lcd$x
    d = ncol(y)
    n = nrow(y)
    logMLE  = lcd$logMLE
    ntriang = nrow(triang)
    slcd <- rep(0,n)

    wt_com = simplex_unit_05_2_nd_pt (d)$coef
    unit_pt_com = simplex_unit_05_2_nd_pt (d)$pt
    wt_com = rep(wt_com,ntriang)
    wt_lcd    = rep(lcd$detA,each=nrow(unit_pt_com))
        
    d_lcd = NULL
    pt    = NULL
    for (j in 1:ntriang) {
        d_lcd = c(d_lcd,exp(crossprod(t(unit_pt_com),logMLE[triang[j,]])))
        pt = rbind(pt,unit_pt_com %*% x[triang[j,],])
    }
    for (j in 1:n) {
        slcd[j] <- sum(wt_com * wt_lcd * d_lcd * dmvnorm(pt,y[j,],A))
    }
    return (slcd)
}

'dslcd_2D' <- function (y, lcd, A) 
{
#  Evaluate the density by Gaussian quadrature
   d = ncol(y)
   n = nrow(y)
   if (d == 2) {
       triang <- lcd$triang
       x <- lcd$x
       logMLE <- lcd$logMLE
       nrows <- nrow(triang)
       slcd <- rep(0,n)
       for(i in 1:n) {            
           for (j in 1:nrows) {
               invA = ginv(A)
               A1 = x[triang[j,1],1] - y[i,1]
               A2 = x[triang[j,1],2] - y[i,2]
               B1 = x[triang[j,2],1] - x[triang[j,1],1] 
               B2 = x[triang[j,2],2] - x[triang[j,1],2] 
               C1 = x[triang[j,3],1] - x[triang[j,1],1] 
               C2 = x[triang[j,3],2] - x[triang[j,1],2] 
               fA = logMLE[triang[j,1]]
               fB = logMLE[triang[j,2]]
               fC = logMLE[triang[j,3]]
               aa = fA - 0.5* A1^2 * invA[1,1] -  0.5 * A2^2 * invA[2,2] - A1 * A2 * invA[1,2]
               bb = -fA + fB - invA[1,1] * B1*A1 - invA[2,2] * B2*A2 -  invA[1,2] * (A1*B2 + B1*A2)
               cc = -fA + fC - invA[1,1] * C1*A1 - invA[2,2] * C2*A2 -  invA[1,2] * (A1*C2 + C1*A2)
               dd = 0.5 * invA[1,1] * B1^2 + 0.5 * invA[2,2] * B2^2 +  invA[1,2] * B1 * B2 
               ee = 0.5 * invA[1,1] * C1^2 + 0.5 * invA[2,2] * C2^2 +  invA[1,2] * C1 * C2 
               ff = - invA[1,2] * (B1*C2 + B2*C1) - invA[1,1] * C1 * B1 - invA[2,2] * C2 * B2
               
               # f is the 1D integration of exp(aa + bb x + cc y - dd x2 - eey2 + ffxy)
               f<-function(x){
                   part1<-exp((cc^2+4*aa*ee+ 2*cc*ff*x + x*(4*bb*ee-4*dd*ee*x+ff^2*x))/4/ee)
                   part2<-pnorm((cc+ff*x)/(ee^0.5)/sqrt(2))*2-1
                   part3<-pnorm((cc+ff*x+2*ee*(x-1))/(ee^0.5)/sqrt(2))*2-1
                   answer = part1 * pi^0.5 * (part2-part3)/2/(ee^0.5)                      
               }
               slcd[i] <- slcd[i] + lcd$detA[j] / 2 / pi /det(A)^0.5 * integrate(f,0,1)$value
           }
       }
       return(slcd)
    }
}

'dslcd_1D' <- function (y, lcd, A) 
{
#  Evaluate the density in one dimension via explicit expression
   d = ncol(y)
   n = nrow(y)
   if (d == 1) {
       triang = lcd$triang
       x = lcd$x
       A = A[1]
       logMLE = lcd$logMLE
       nrows = nrow(triang)
       slcd = rep(0,n)
       for(i in 1:n) {            
           for (j in 1:nrows) {
               x1 = x[triang[j,1]]
               x2 = x[triang[j,2]]
               phi1 = logMLE[triang[j,1]]
               phi2 = logMLE[triang[j,2]]
               s = (phi2-phi1)/(x2-x1)    
               slcd[i] = slcd[i] + abs(exp(phi1+s*(A*s/2-x1+y[i]))*(pnorm(x2,y[i]+A*s,sqrt(A))-pnorm(x1,y[i]+A*s,sqrt(A))))
            }
       }
       return(slcd)
    }
}
