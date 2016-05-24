expertprobs <-
function(fit, x, d = "best", ex = 1){
	
	if(d == "best"){
		ssq <- fit$ssq[ex, 1 - is.na(fit$ssq[ex,])]
		best.index <- which(ssq == min(ssq))[1]
	}
	index<-switch(which(d==c("normal", "t", "gamma", "lognormal", "logt","beta", "best")), 1, 2, 3, 4, 5, 6, best.index)
	
		
	if(index==1){
		px <- pnorm(x, fit$Normal[ex,1], fit$Normal[ex,2]) 		
	}
	
	if(index==2){
		px <- pt((x - fit$Student.t[ex,1])/fit$Student.t[ex,2], fit$Student.t[ex,3])	
	}
	
	if(index==3){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		px <- pgamma(x - xl, fit$Gamma[ex,1], fit$Gamma[ex,2])  
	}
	
	if(index==4){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		px <- plnorm(x - xl, fit$Log.normal[ex,1], fit$Log.normal[ex,2]) 
	}	
	
	if(index==5){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		px <- pt( (log(x - xl) - fit$Log.Student.t[ex,1]) / fit$Log.Student.t[ex,2], fit$LogStudent.t[ex,3]) 
	}
		
	if(index==6){
		xl <- fit$limits[ex,1]
		xu <- fit$limits[ex,2]
		if(xl == -Inf){xl <- 0}
		if(xu == Inf){xu <- 1}
		px <-  pbeta( (x - xl) / (xu - xl), fit$Beta[ex,1], fit$Beta[ex,2])
	}
	
  
  
px	
	
}
