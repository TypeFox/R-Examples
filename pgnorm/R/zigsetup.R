zigsetup <-
function(p,n,tol){

# A function implemented by Steve Kalke

# Description: 
# Set's up the ziggurat for the one-sided, p-generalized normal distribution
# by returning a vector with the rightmost endpoints of the first n-1 rectangles


# Arguments: 
# p-   a positiv constant (default: p=2)
# n-   a positive integer, expressing the number of rectangles(default: n=2^8)
# tol- a positive constant, defining the accuracy of evaluation (default: tol=10^(-9))

if(missing(p)){p<-2}

if(p<0){stop("p has to be positive")}

if(missing(n)){n<-2^8}

if(round(n)!=n |n<=0){stop("invalid argument for positive integer n")}

if(missing(tol)){tol<-10^(-9)}


x<-rep(1,n)
            # evaluating the rightmost endpoint for the search-interval of x[n] 
# if the area of a rectangle is less than 1/n, x1 exceeds x[n]
x1<-10  
while( x1*2*dpgnorm(x1,p)+1-pgamma(((x1)^p)/p,1/p)>= 1/n ){x1<-10*x1} #increase x1, until the area of a rectangle is greater than 1/n

x0<-0#leftmost endpoint for the search-interval of x[n] 

while(abs(x[1])>tol){
x[n]<-(x0+x1)/2  # nested intervals method
x[1]<-(-1)
v<-x[n]*2*dpgnorm(x[n],p)+1-pgamma((x[n]^p)/p,1/p) #area of the rectangles defined by n and x[n]
for (i in seq(n-1,1,-1)){
if( gamma(1/p)*(p^(1/p-1))*(v/x[i+1]+2*dpgnorm(x[i+1],p)) > 1)break #argument out of domain of the inverse function
#-> x[n] too small
if(i==1){x[1]<-x[2]-v/2/(dpgnorm(0,p)-dpgnorm(x[2],p))}
else{x[i]<-(-p*log(gamma(1/p)*(p^(1/p-1))*(v/x[i+1]+2*dpgnorm(x[i+1],p))))^(1/p) }
}
if(x[1]>tol){x1<-x[n]}# x[n] too big
if(x[1]<(-1*tol)){x0<-x[n]}# x[n] too small

}
return(x[2:n])
}
