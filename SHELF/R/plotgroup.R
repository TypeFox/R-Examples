plotgroup <-
function(fit, pl, pu, d = "best", bw = T){
	
  
  
	n.experts <- nrow(fit$vals)
	x <- matrix(0, 200, n.experts)
	fx <- x
	
	if(bw == T){colvec <- rep(1, n.experts)}else{colvec <- c(1:n.experts)}
	
	for(i in 1:n.experts){
		densitydata <- expertdensity(fit, d, ex = i, pl, pu)
		x[,i] <- densitydata$x
		fx[,i] <-densitydata$fx
	}
	
	par(ps=15)
	par(mar = c(5.1, 5.1, 4.1, 2.1))
	
	plot(x[,1], fx[,1], type = "l", lwd = 2, xlab = "x", ylab = expression(f[X](x)), xlim = c(min(x), max(x)), ylim = c(min(fx), max(fx)), col=colvec[1])
	for(i in 2:n.experts){
		lines(x[,i], fx[,i], lwd = 2, lty = i, col = colvec[i])
	}
	
	xmax <- x[fx==max(fx)]
	
	if(pl == -Inf){pl <- qnorm(0.001, fit$Normal[1,1], fit$Normal[1,2])}
	if(pu == Inf){pu <- qnorm(0.999, fit$Normal[1,1], fit$Normal[1,2])}
  
	if(xmax[1] < mean (c(pl,pu))){
		legend("topright", legend=c(1:n.experts), lty=c(1:n.experts), lwd=rep(2,n.experts), col = colvec)}else{
		legend("topleft", legend=c(1:n.experts), lty=c(1:n.experts), lwd=rep(2,n.experts), col = colvec)
	}
		
}
