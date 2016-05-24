expertquantiles <-
function(fit, q, d = "best", ex = 1){
	
	if(d == "best"){
		ssq <- fit$ssq[ex, 1 - is.na(fit$ssq[ex,])]
		best.index <- which(ssq == min(ssq))[1]
	}
	index<-switch(which(d==c("normal", "t", "gamma", "lognormal", "logt","beta", "best")), 1, 2, 3, 4, 5, 6, best.index)
	
		
	if(index==1){
		qx <- qnorm(q, fit$Normal[ex,1], fit$Normal[ex,2]) 		
	}
	
	if(index==2){
		qx <- fit$Student.t[ex,1] + fit$Student.t[ex,2] * qt(q, fit$Student.t[ex,3])	
	}
	
	if(index==3){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		qx <- xl + qgamma(q, fit$Gamma[ex,1], fit$Gamma[ex,2])  
	}
	
	if(index==4){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		qx <- xl + qlnorm(q, fit$Log.normal[ex,1], fit$Log.normal[ex,2]) 
	}	
	
	if(index==5){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		qx <- xl + exp(fit$Log.Student.t[ex,1] + fit$Log.Student.t[ex,2] * qt( q , fit$LogStudent.t[ex,3])) 
	}
		
	if(index==6){
		xl <- fit$limits[ex,1]
		xu <- fit$limits[ex,2]
		if(xl == -Inf){xl <- 0}
		if(xu == Inf){xu <- 1}
		qx <- xl + (xu - xl) * qbeta(q, fit$Beta[ex,1], fit$Beta[ex,2])
	}
	
qx	
	
}
