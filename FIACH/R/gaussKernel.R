gaussKernel<-function(fwhm=8,ksize=5,voxsize=1){
sigma <- fwhm^2 / (8 * log(2))
if(ksize%%2==0){ksize<-ksize+1}
surround<-floor(ksize/2)
x<- (-surround:surround)*voxsize
kern<-dnorm(x = x,mean = 0,sd = sqrt(sigma))
kern<-kern/sum(kern)
return(kern)
}

