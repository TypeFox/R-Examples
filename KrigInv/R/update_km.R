update_km <- function(model, NewX,NewY,NewX_AllreadyExist=FALSE,CovReEstimate=FALSE,new.noise.var=NULL,kmcontrol=NULL,F.newdata=NULL){

	#model: a km object
	#NewX: a matrix containing the new points of experiments
	#NewY: a matrix containing the function values on the points NewX
	#NewX_AllreadyExist: a boolean: do the NewX points allready exist in the design of experiment model@X ?
	#CovReEstimate: a boolean: do you want to re-estimate the covariance parameters with the new observations ?
	#new.noise.var: a vector containing the noise variance at each new observations
	#F.newdata: an optional matrix containing the value of the trend at the new locations (to avoid a
		#VERY expensive "model.matrix" call)
	
	#This function updates a km model in 3 possible different ways:
		# - with new points and new responses WITH covariance parameters re-estimation
		# - with new points and new responses WITH NO covariance parameters re-estimation
		# - with EXISTING points and modified responses for these points, with NO covariance parameters re-estimation
		
	if(NewX_AllreadyExist == FALSE){
		#NewX <- matrix(NewX,ncol=model@d)
		
		model@X <- rbind(model@X, as.matrix(NewX))
		model@y <- rbind(model@y, as.matrix(NewY))
		
		#Consistency tests between model@noise.var and new.noise.var
		if ((length(model@noise.var) != 0) && (is.null(new.noise.var))) {
			#noisy initial observations and no noise for new observations
			if (is.null(nrow(NewX))) {TheNew.noise.var=0}
			else {TheNew.noise.var <- rep(0,nrow(NewX))}								
			model@noise.var <- c(model@noise.var,TheNew.noise.var)
		}
		if ((length(model@noise.var) != 0) && (!is.null(new.noise.var))){
			#noisy initial observations and noisy new observations
			model@noise.var <- c(model@noise.var,new.noise.var)
		}
		if ((length(model@noise.var) == 0) && (!is.null(new.noise.var))){
			#noise free initial observations and noisy new observations
			if (any(new.noise.var!=0)){	#when previous noise are 0 AND all the new noises are 0 we are still noise free
				noise.var.init <- rep(0,model@n)
				model@noise.var <- c(noise.var.init,new.noise.var)
			}
		}
		
		if(CovReEstimate){

			#case 1: new points, covariance parameter re-estimation (provided model@param.estim == true)
			if (model@param.estim) {
				#case 1a: here we re-estimate the covariance parameter
				#default values for the cov param estimation parameters - when they are not provided
				if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
				if (length(model@penalty==0)) kmcontrol$penalty <- NULL 
				if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
				if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- covparam2vect(model@covariance)
				if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
				
				model <- km(formula=model@trend.formula, design=model@X,response=model@y,
							covtype=model@covariance@name, lower=model@lower,upper=model@upper,
							noise.var=model@noise.var,penalty=kmcontrol$penalty, optim.method=kmcontrol$optim.method,
							parinit=kmcontrol$parinit, control=kmcontrol$control, gr=model@gr)
			}
			else{
				#case 1b: no re-estimation. Keeping the current covariance parameters
				
				#solution 1 (in comments because it's slower)
				#model <- km(formula=model@trend.formula, design=model@X,response=model@y,
				#			covtype=model@covariance@name,coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),coef.var=model@covariance@sd2, 
				#			noise.var=model@noise.var)
					
				#solution 2	(faster)
				model@n <- nrow(model@X)
				#model@F <- model.matrix(model@trend.formula, data=data.frame(model@X))
				if (is.null(F.newdata)) {model@F<-trendMatrix.update(model,Xnew=data.frame(NewX))
				} else model@F<-rbind(model@F,F.newdata)
				
				model <- computeAuxVariables_update(model)
			
			}
		}
		else{
			#case 2: new points, no covariance parameter re-estimation. keeping the current covariance parameters
			
			#solution 1 (in comments because it's slower)
			#model <- km(formula=model@trend.formula, design=model@X,response=model@y,
			#			covtype=model@covariance@name,coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),coef.var=model@covariance@sd2, 
			#			noise.var=model@noise.var)
						
			#solution 2	(faster)
			
			model@n <- nrow(model@X)
			#model@F <- model.matrix(model@trend.formula, data=data.frame(model@X))
			if (is.null(F.newdata)) {model@F<-trendMatrix.update(model,Xnew=data.frame(NewX))
			} else model@F<-rbind(model@F,F.newdata)
			
			model <- computeAuxVariables_update(model)
			
		}
	}
	else{
		#case 3: existing points with a modified response, 
		#No need to recalculate the Cholesky but the other variables are still needed (z, M) for prediction. 
		#No covariance param re-estimation in this case
		#we assume that the k new values given modify the k last points of the design model@X
		
		if (is.null(nrow(NewX))) numrow=1
		else numrow=nrow(NewX)
		for(i in 1:numrow) model@y[model@n-numrow+i] <- NewY[i]
		model <- computeAuxVariables_noChol(model)
	}
	return(model)
}
