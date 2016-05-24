dgam <-
function(y,mu,delta){
    dgamma(y,shape=1/delta,rate=1/(delta*mu))  
}
