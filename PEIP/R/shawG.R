###  Parameter Estimation and Inverse Problems, 2nd edition, 2011
###  by R. Aster, B. Borchers, C. Thurber
### Generate the m by n system matrix for the Shaw problem 
###  
###  G = shaw.G(n,m) 
###  
###  
###  m and n must be even. 
 
###  Reference: C. B. Shaw, Jr., "Improvements of the resolution of 
###  an instrument by numerical solution of an integral equation", 
###  J. Math. Anal. Appl. 37m 83-112, 1972. 
shawG <-function(m,n)  {  
 
if( (m %% 2) !=0) { print('m must be even'); return()  } # 
if( (n%%2)!=0) {  print('n must be even'); return()     } # 
 
###  Initialization. 
delta = pi/n
G = matrix(rep(0,times=m*n) ,m,n)

vm=((1:m)-1/2)*pi/m-pi/2
vn=((1:n)-1/2)*pi/n-pi/2

for (  i in 1:m ) {
    for (  j in 1:n ) {
        s=vm[i]
        th=vn[j]
        x=pi*(sin(s)+sin(th))
        if(x==0)
          {
            G[i,j]=(cos(s)+cos(th))^2
          }
        else
          {
        G[i,j]=(cos(s)+cos(th))^2*(sin(x)/x)^2
      }
           } 
  }



G=G * delta


return(G)
}
