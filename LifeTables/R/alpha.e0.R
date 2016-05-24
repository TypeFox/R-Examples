alpha.e0 <-
function(pattern, e0.target, sex="female"){
	
	
		ex.alpha <- function(alpha){
		
		pred <- mortmod(pattern=pattern, alpha=alpha, sex=sex)
		e0 <- lt.mx(nmx=exp(pred),age=c(0,1,seq(5,110,5)), sex=sex)$e0 
		ex.diff <- abs(e0.target-e0)
		return(ex.diff)
	}
	optimize(f=ex.alpha, interval=c(-5,3), tol=.01)$minimum
	}

