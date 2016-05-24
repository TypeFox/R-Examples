pgam <-
function(y,mu,delta){
    pgamma(y,shape=1/delta,rate=1/(delta*mu))  
}
