`waller` <-
function(K,q,f,Fc) {
a0<-1
b0<-20
for (i in 1:50) {
t<-(b0+a0)/2
g<-function(x,q,f,Fc) x^(-(q+3)/2)*(f+q*Fc/x)^(-(f+q-1)/2)
b<-function(x,q,f,Fc) sqrt((f+q)*(x-1)/(f*x+q*Fc) )
h<-function(z,q,f) ((f+q+z^2)/(f+q-1))*dt(z,q+f)+z*pt(z,q+f)
n0 <- function(x) sqrt(x-1)*g(x,q,f,Fc)*h( t*b(x,q,f,Fc),q,f )
d0<- function(x) sqrt(x-1)*g(x,q,f,Fc)*h( -t*b(x,q,f,Fc),q,f )
n1 <- function(x) sqrt(x-1)*g(x,q,f,Fc)*h( a0*b(x,q,f,Fc),q,f )
d1<- function(x) sqrt(x-1)*g(x,q,f,Fc)*h( -a0*b(x,q,f,Fc),q,f )

# Teoria k = int(numerador,1,inf)/int(denominador,1 inf)
#------------------------
IN0<-integrate(n0,1,Inf)$value
ID0<-integrate(d0, 1, Inf)$value
IN1<-integrate(n1,1,Inf)$value
ID1<-integrate(d1, 1, Inf)$value

if( (K-IN0/ID0)*(K-IN1/ID1) <= 0) b0<- t
if( (K-IN0/ID0)*(K-IN1/ID1) > 0 ) a0<- t
if ( abs(b0-a0) <= 5.0e-4 ) break
}
return(round(t,3))
}

