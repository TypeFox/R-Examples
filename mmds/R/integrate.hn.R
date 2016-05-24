integrate.hn<-function(keysc,upper){
   # integarte the half-Normal detection function...

   sqrt(2*pi)*keysc*(pnorm(upper,sd=keysc)-pnorm(0,sd=keysc))

}
