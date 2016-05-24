.logMCMC.post <- function (x,n,theta){

 if(length(x)!=length(n) | length(x)!=length(theta)){
	print("data lengths do not match!")
	stop()
 }
 sum(dbinom(x,n,prob=theta,log=T)/length(x))

}
