rgam <-
function(n,mu,delta){
    rgamma(n,shape=1/delta,rate=1/(delta*mu))  
}
