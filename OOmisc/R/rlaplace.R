rlaplace <-
function(n,beta){
if(beta<=0) stop("Beta should be positive") 
u=runif(n)
x=ifelse(u<0.5,beta*log(2*u),-beta*log(2*(1-u)))
if(beta>0) x
}

