## Growth curve modeling
## May 23, 2012
## Johnny Zhang, Jack McArdle, & Aki Hamagami

ramLCM<-function(data, outcome, model=c('all','no','linear','quadratic','latent'), basis=0:(length(outcome)-1), predictor, equal.var=TRUE, digits=3, ram.out=FALSE, ...){
	if (missing(data)) stop("No data was provided!")
	if (missing(outcome)) stop("Missing outcome variables!")
	if (!is.data.frame(data)) stop("The provided data set is not a data frame!")
	if (length(model)>1) model="all"
	
	varname<-names(data)
	n<-length(outcome)
	if (!missing(predictor)) np<-length(predictor)
	## general lavaan input for all models
	
	## 1. the no growth model
	model.no<-paste("level =~ 1*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+1*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	## add predictors
	if (!missing(predictor)){
		model.no<-paste(model.no, "level~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[i]])
		}
		model.no<-paste(model.no, "\n")
	}
	if (equal.var){
		for (i in 1:n){
			model.no<-paste(model.no, varname[outcome[i]], "~~(vare)*", varname[outcome[i]], "\n")
		}
	}
	growth.no<-model.no
	
	## 2. the linear growth model
	model.no<-NULL
	model.no<-paste("level =~ 1*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+1*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	model.no<-paste(model.no, "slope =~ ",basis[1], "*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+", basis[i], "*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	## add predictors
	if (!missing(predictor)){
		model.no<-paste(model.no, "level~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
		
		model.no<-paste(model.no, "slope~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
	}
	if (equal.var){
		for (i in 1:n){
			model.no<-paste(model.no, varname[outcome[i]], "~~(vare)*", varname[outcome[i]], "\n")
		}
	}
	growth.linear<-model.no
	
	## 3.  the latent growth model
	model.no<-NULL
	model.no<-paste("level =~ 1*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+1*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	model.no<-paste(model.no, "slope =~ ",basis[1], "*", varname[outcome[1]])
	
	for (i in 2:(n-1)){
		model.no<- paste(model.no, "+start(", basis[i], ")*", varname[outcome[i]])
	}
	model.no<- paste(model.no, "+", basis[n], "*", varname[n])
	
	model.no<-paste(model.no, "\n")
	
	## add predictors
	if (!missing(predictor)){
		model.no<-paste(model.no, "level~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
		
		model.no<-paste(model.no, "slope~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
	}
	if (equal.var){
		for (i in 1:n){
			model.no<-paste(model.no, varname[outcome[i]], "~~(vare)*", varname[outcome[i]], "\n")
		}
	}
	growth.latent<-model.no
	
	## 4.  the quadratic growth model
	model.no<-NULL
	model.no<-paste("level =~ 1*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+1*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	model.no<-paste(model.no, "slope =~ ",basis[1], "*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+", basis[i], "*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	model.no<-paste(model.no, "quadratic =~ ", basis[1]^2, "*", varname[outcome[1]])
	
	for (i in 2:n){
		model.no<- paste(model.no, "+", basis[i]^2, "*", varname[outcome[i]])
	}
	
	model.no<-paste(model.no, "\n")
	
	## add predictors
	if (!missing(predictor)){
		model.no<-paste(model.no, "level~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
		
		model.no<-paste(model.no, "slope~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
		
		model.no<-paste(model.no, "quadratic~", varname[predictor[1]])
		if (np>1){
			for (j in 2:np) model.no<- paste(model.no, "+", varname[predictor[j]])
		}
		model.no<-paste(model.no, "\n")
	}
	if (equal.var){
		for (i in 1:n){
			model.no<-paste(model.no, varname[outcome[i]], "~~(vare)*", varname[outcome[i]], "\n")
		}
	}
	growth.quadratic<-model.no

	## fit the model using lavaan
	if (model=="no"){
		results.no<-growth(model=growth.no, data=data,...)
		ram.no<-lavaan2ram(results.no, ram.out=ram.out)
		summary(results.no, fit.measures=TRUE)
		invisible((list(lavaan=list(no=results.no), ram=list(no=ram.no), model=list(no=growth.no, linear=growth.linear, latent=growth.latent, quadratic=growth.quadratic), data=data, type=model)))
	}else if (model=="linear"){
		results.linear<-growth(model=growth.linear, data=data,...)
		ram.linear<-lavaan2ram(results.linear, ram.out=ram.out)
		summary(results.linear,fit.measures=TRUE)
		invisible((list(lavaan=list(linear=results.linear), ram=list(linear=ram.linear), model=list(no=growth.no, linear=growth.linear, latent=growth.latent, quadratic=growth.quadratic), data=data, type=model)))
	}else if (model=="latent"){
		results.latent<-growth(model=growth.latent, data=data,...)
		ram.latent<-lavaan2ram(results.latent, ram.out=ram.out)
		summary(results.latent,fit.measures=TRUE)
		invisible((list(lavaan=list(latent=results.latent), ram=list(latent=ram.latent), model=list(no=growth.no, linear=growth.linear, latent=growth.latent, quadratic=growth.quadratic), data=data, type=model)))
	}else if (model=="quadratic"){
		results.quadratic<-growth(model=growth.quadratic, data=data,...)
		ram.quadratic<-lavaan2ram(results.quadratic, ram.out=ram.out)
		summary(results.quadratic,fit.measures=TRUE)
		invisible((list(lavaan=list(quadratic=results.quadratic), ram=list(quadratic=ram.quadratic), model=list(no=growth.no, linear=growth.linear, latent=growth.latent, quadratic=growth.quadratic), data=data, type=model)))
	}else if (model=="all"){
		
		results.no<-growth(model=growth.no, data=data,...)		
		cat("\n\n======================================\nResults from the no growth curve model\n======================================\n\n")
		ram.no<-lavaan2ram(results.no, ram.out=ram.out)
		summary(results.no,fit.measures=TRUE)
		
		cat("\n\n==========================================\nResults from the linear growth curve model\n==========================================\n\n")
		results.linear<-growth(model=growth.linear, data=data,...)
		ram.linear<-lavaan2ram(results.linear, ram.out=ram.out)
		summary(results.linear,fit.measures=TRUE)
		
		cat("\n\n==========================================\nResults from the latent growth curve model\n==========================================\n\n")
		results.latent<-growth(model=growth.latent, data=data,...)
		ram.latent<-lavaan2ram(results.latent, ram.out=ram.out)
		summary(results.latent,fit.measures=TRUE)
		
		cat("\n\n=============================================\nResults from the quadratic growth curve model\n=============================================\n\n")
		results.quadratic<-growth(model=growth.quadratic, data=data,...)
		ram.quadratic<-lavaan2ram(results.quadratic, ram.out=ram.out)
		summary(results.quadratic,fit.measures=TRUE)
		
		## Compare models
		fit.no<-fitMeasures(results.no)
		fit.linear<-fitMeasures(results.linear)
		fit.latent<-fitMeasures(results.latent)
		fit.quadratic<-fitMeasures(results.quadratic)
		cat("\n\n====================================================\nFit Statistics and Fit Indices for Model Comparisons\n====================================================\n\n")
		all.fit<-cbind(fit.no,fit.linear, fit.latent, fit.quadratic)
		chisq.per.df<-all.fit[1,]/all.fit[2,]
		n.all<-nrow(all.fit)
		all.fit<-rbind(all.fit[1:2, ], chisq.per.df, all.fit[3:n.all, ])		
		colnames(all.fit)<-c("No","Linear","Latent","Quadratic")
		print(all.fit,digits=digits)
		invisible((list(lavaan=list(no=results.no, linear=results.linear, latent=results.latent, quadratic=results.quadratic), ram=list(no=ram.no, linear=ram.linear, latent=ram.latent, quadratic=ram.quadratic), fit=all.fit, model=list(no=growth.no, linear=growth.linear, latent=growth.latent, quadratic=growth.quadratic), data=data, type=model)))
	}else{
		stop("Wrong model option was used!")
	}
}

