plotlinearpool <-
function(fit, xl, xu, ql = NA, qu = NA, d = "best", ind = T, w = 1){
	
	n.experts <- nrow(fit$vals)
	x <- matrix(0, 200, n.experts)
	fx <- x
  if(min(w)<0 | max(w)<=0){stop("expert weights must be non-negative, and at least one weight must be greater than 0.")}
  
	if(length(w)==1){
	  w <- rep(w, n.experts)
	}
  
	weight <- matrix(w/sum(w), 200, n.experts, byrow = T)
 
	for(i in 1:n.experts){
		densitydata <- expertdensity(fit, d, ex = i, xl, xu)
		x[,i] <- densitydata$x
		fx[,i] <-densitydata$fx 
	}
	
	fx.lp <- apply(fx*weight,1,sum)
	
	par(ps=15)
	par(mar = c(5.1, 5.1, 4.1, 2.1))
	
	plot(x[,1], fx.lp, type = "l", lwd = 4, xlab = "x", ylab = expression(f[X](x)), xlim = c(min(x), max(x)), ylim = c(min(fx), max(fx)), col="red")
	
	if(ind == T){
	
		for(i in 1:n.experts){
			lines(x[,i], fx[,i], lwd = 2, lty = 2)
		}
	
		xmax <- x[fx==max(fx)]
	
		if(xmax[1] < mean (c(xl,xu))){
			legend("topright", legend=c("individual", "pooled"), lty=c(2,1), lwd=c(2,4), col = c(1,2))}else{
			legend("topleft", legend=c("individual", "pooled"), lty=c(2,1), lwd=c(2,4), col = c(1,2))	}
	}
		
	if(is.na(ql) == F){
		x.q1 <- qlinearpool(fit, ql, d, w)
		if(x.q1 > xl){
		  x1 <-seq(from = xl, to = x.q1 , length = ceiling((x.q1-xl)/(xu - xl)*200))
		  lines(x1, approx(x = x[,1], y = fx.lp, xout = x1)$y, type="h", col = "red" )
		}
	}
	
	if(is.na(qu) == F){
		x.q2 <- qlinearpool(fit, qu, d, w)
		if(x.q2 < xu){
		  x2<-seq(from = x.q2, to = xu, length = ceiling((xu - x.q2)/(xu - xl)*200))
		  lines(x2, approx(x = x[,1], y = fx.lp, xout = x2)$y , type="h", col = "red" )
	  }
	}
			
}
