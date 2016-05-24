icproportion<-function(prop=prop, n=n, alpha=0.05){
   inter=alpha/2#seuil de confiance
   z=qnorm(1-inter)
   basic=(2*n*prop+z^2)/(2*(n+z^2))
   complement=(z*sqrt(z^2+4*n*prop*(1-prop)))/(2*(n+z^2))
   lower =basic-complement
   upper =prop + complement
   return(c(lower, upper))
}