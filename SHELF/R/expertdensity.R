expertdensity <-
function(fit, d = "best", ex = 1, pl, pu, nx = 200){
	
  if(pl == -Inf){pl <- qnorm(0.001, fit$Normal[ex,1], fit$Normal[ex,2])}
  if(pu == Inf){pu <- qnorm(0.999, fit$Normal[ex,1], fit$Normal[ex,2])}
  
	x <- seq(from = pl, to = pu, length = nx)
	
	if(d == "best"){
		ssq <- fit$ssq[ex, 1 - is.na(fit$ssq[ex,])]
		best.index <- which(ssq == min(ssq))[1]
	}
	index<-switch(which(d==c("normal", "t", "gamma", "lognormal", "logt","beta", "hist", "best")), 1, 2, 3, 4, 5, 6, 7, best.index)
	
		
	if(index==1){
		fx <- dnorm(x, fit$Normal[ex,1], fit$Normal[ex,2]) 		
	}
	
	if(index==2){
		fx <- dt((x - fit$Student.t[ex,1])/fit$Student.t[ex,2], fit$Student.t[ex,3])/fit$Student.t[ex,2]
	}
	
	if(index==3){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		fx <- dgamma(x - xl, fit$Gamma[ex,1], fit$Gamma[ex,2])  
	}
	
	if(index==4){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		fx <- dlnorm(x - xl, fit$Log.normal[ex,1], fit$Log.normal[ex,2]) 
	}	
	
	if(index==5){
		xl <- fit$limits[ex,1]
		if(xl == -Inf){xl <- 0}
		fx <- dt( (log(x - xl) - fit$Log.Student.t[ex,1]) / fit$Log.Student.t[ex,2], fit$Log.Student.t[ex,3]) / ((x - xl) * fit$Log.Student.t[ex,2])
    fx[is.nan(fx)]<-0
    
	}
		
	if(index==6){
		xl <- fit$limits[ex,1]
		xu <- fit$limits[ex,2]
		if(xl == -Inf){xl <- 0}
		if(xu == Inf){xu <- 1}
		fx <-  1/(xu - xl) * dbeta( (x - xl) / (xu - xl), fit$Beta[ex,1], fit$Beta[ex,2])
	}

	if(index==7){

	  fx <- dhist(x,fit$vals[ex,], fit$probs[ex,])}

 
list(x = x, fx = fx)	
	
}
