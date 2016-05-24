riskscoreci <-
function(x1,n1,x2,n2,conf.level)
{
  z =  abs(qnorm((1-conf.level)/2))
  if ((x2==0) &&(x1==0)){
    ul = Inf
    ll = 0
    }
  else{
     a1 =  n2*(n2*(n2+n1)*x1+n1*(n2+x1)*(z^2))
     a2 = -n2*(n2*n1*(x2+x1)+2*(n2+n1)*x2*x1+n1*(n2+x2+2*x1)*(z^2))
     a3 = 2*n2*n1*x2*(x2+x1)+(n2+n1)*(x2^2)*x1+n2*n1*(x2+x1)*(z^2)
     a4 = -n1*(x2^2)*(x2+x1)
     b1 = a2/a1
     b2 = a3/a1
     b3 = a4/a1
     c1 = b2-(b1^2)/3
     c2 = b3-b1*b2/3+2*(b1^3)/27
     ceta = acos(sqrt(27)*c2/(2*c1*sqrt(-c1)))
     t1 = -2*sqrt(-c1/3)*cos(pi/3-ceta/3)
     t2 = -2*sqrt(-c1/3)*cos(pi/3+ceta/3)
     t3 = 2*sqrt(-c1/3)*cos(ceta/3)
     p01 = t1-b1/3
     p02 = t2-b1/3
     p03 = t3-b1/3
     p0sum = p01+p02+p03
     p0up = min(p01,p02,p03)
     p0low = p0sum-p0up-max(p01,p02,p03)

     if( (x2==0) && (x1!=0) ){
        ll = (1-(n1-x1)*(1-p0low)/(x2+n1-(n2+n1)*p0low))/p0low
        ul = Inf
       }
     else if( (x2!=n2) && (x1==0)){
        ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up
        ll = 0
        }
     else if( (x2==n2) && (x1==n1)){
         ul = (n2+z^2)/n2
         ll =  n1/(n1+z^2)
        }
     else if( (x1==n1) || (x2==n2) ){
         if((x2==n2) && (x1==0)) { ll = 0 }
         if((x2==n2) && (x1!=0)) {
           phat1  = x2/n2
           phat2  =  x1/n1
           phihat = phat2/phat1
           phil = 0.95*phihat
           chi2 = 0
           while (chi2 <= z){
             a = (n2+n1)*phil
             b = -((x2+n1)*phil+x1+n2)
             c = x2+x1
             p1hat = (-b-sqrt(b^2-4*a*c))/(2*a)
             p2hat = p1hat*phil
             q2hat = 1-p2hat
             var = (n2*n1*p2hat)/(n1*(phil-p2hat)+n2*q2hat)
             chi2 = ((x1-n1*p2hat)/q2hat)/sqrt(var)
             ll = phil
             phil = ll/1.0001}}
         i = x2
         j = x1
         ni = n2
         nj = n1
         if( x1==n1 ){
            i = x1
            j = x2
            ni = n1
            nj = n2
         }
         phat1  = i/ni
         phat2  =  j/nj
         phihat = phat2/phat1
         phiu = 1.1*phihat
         if((x2==n2) && (x1==0)) {
            if(n2<100) {phiu = .01}
            else {phiu=0.001}
           }
         chi1 = 0
         while (chi1 >= -z){
         a = (ni+nj)*phiu
         b = -((i+nj)*phiu+j+ni)
         c = i+j
         p1hat = (-b-sqrt(b^2-4*a*c))/(2*a)
         p2hat = p1hat*phiu
         q2hat = 1-p2hat
         var = (ni*nj*p2hat)/(nj*(phiu-p2hat)+ni*q2hat)
         chi1  = ((j-nj*p2hat)/q2hat)/sqrt(var)
         phiu1 = phiu
         phiu = 1.0001*phiu1
         }

         if(x1==n1) {
          ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up
          ll = 1/phiu1
         }
         else{ ul = phiu1}
       }

     else{
     ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up
     ll = (1-(n1-x1)*(1-p0low)/(x2+n1-(n2+n1)*p0low))/p0low
      }
   }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint)
  class(rval) <- "htest"
  return(rval)
}

