`norm.gauss.hermite` <-
function(n){
  rule <- gaussHermiteData(n)
  r <- cbind(rule$x,rule$w)  
  r[,1] <- r[,1]*sqrt(2)
	r[,2] <- r[,2]/sqrt(pi)
#	r[,2] <- r[,2]/sum(r[,2])
	return(r)
}
