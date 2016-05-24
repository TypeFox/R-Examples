# x=(mu1, delta, sigma1, sigma2); delta=mu1-mu2
LoglikNormalDeltaCens <- function(x, data1, lowerbound1, vdelta1, data2, lowerbound2, vdelta2){
	if(x[3] > 0 & x[4] > 0){
		initial1 <- c(x[1],x[3])
		initial2 <- c(x[1]-x[2],x[4])
		loglik <- LoglikNormalCens(x=initial1, data=data1, lowerbound=lowerbound1, vdelta=vdelta1) + LoglikNormalCens(x=initial2, data=data2, lowerbound=lowerbound2, vdelta=vdelta2)
	}
     
	if(x[3] <= 0 | x[4] <= 0){ loglik <- -Inf }
     if(loglik == Inf){loglik <- -100000}

     return(loglik)
}










#