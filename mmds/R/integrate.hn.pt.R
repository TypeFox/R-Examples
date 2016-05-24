integrate.hn.pt<-function(keysc,upper){
   # integarte the half-Normal detection function for point transects
   2*pi*(keysc^2)*(1-exp(-(upper^2)/(2*keysc^2)))
}
