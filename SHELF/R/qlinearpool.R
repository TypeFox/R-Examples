#' @export
#' 
qlinearpool <-
function(fit, q, d = "best", w = 1){
	
	n.experts <- nrow(fit$vals)
	qx.individual <- matrix(0, length(q), n.experts)
	
	for(i in 1:n.experts){
		qx.individual[,i] <- expertquantiles(fit, q, d, ex = i)
	}
	
	n.q <- length(q)
	qx<-rep(0, n.q)
	
	for(i in 1:n.q){
		x <- seq(from = 0.999*min(qx.individual[i,]), 1.001*max(qx.individual[i,]), length = 100)
		px <- plinearpool(fit, x, d, w)
		qx[i] <- approx(x = px, y = x, xout = q[i])$y 
	}
qx			
}
