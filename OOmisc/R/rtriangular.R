rtriangular <-
function(n,a,b){
if(b<=a) stop("b should be greater than a")  
u=runif(n)
x=ifelse(u<0.5,2*a+(b-a)*sqrt(2*u),2*b-(b-a)*sqrt(2*(1-u)))
if(b>a) x
}

