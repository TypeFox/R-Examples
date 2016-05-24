logitinv=function(x) 1/(1+exp(-x))
logit=function(p) log(p/(1-p))

####################################
############ Creation ##############
####################################

##########  MSM.fitted   ##############

########
#####Class definition
setClass("MSM.fitted",
	representation(           
		CondMean	= "matrix",
		error		= "matrix",
		Likel		= "matrix",
		margLik	= "matrix",
		filtProb	= "matrix",
		smoProb	= "matrix",
		smoTransMat	= "list",
		logLikel	= "numeric"
	)
)

####################################
##### Construction, get and set ####
####################################

##### Get
setMethod(
	f="[",
	signature=c("MSM.fitted","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			CondMean	= return(x@CondMean),
			error		= return(x@error),
			Likel 	= return(x@Likel),
			margLik 	= return(x@margLik),
			filtProb 	= return(x@filtProb),
			smoProb 	= return(x@smoProb),
			states	= return(apply(x@smoProb,1,which.max)),
			smoTransMat	= return(x@smoTransMat),
			logLikel 	= return(x@logLikel),
			stop("Error:",i,"is not a MSM.fitted slot")
		)
	}
)


#############################################
#############################################
#############################################
#############################################
#############################################
#############################################

############  MSM   ################

####################################
############ Creation ##############
####################################

##########
#####Validity function
validMSM=function(object){
	if(object@k<=0){
		stop("A negative k is not possible!")
	}else{}
}

########
#####Class definition

setClass(
	Class="MSM",
	representation=representation(
		call 	= "call",
		k      	= "numeric",
		switch 	= "logical",
		p	= "numeric",
		Fit 	= "MSM.fitted"
	)
)

#############################################
#############################################
#############################################
#############################################
#############################################
#############################################

#########  MSM.linear   ############

####################################
############ Creation ##############
####################################

##########
#####Validity function
validMSM.linear=function(object){
	if(min(object@transMat)<0 | max(object@transMat)>1) {
		stop("A probability must be a value between 0 and 1!")
	}
	if(object@p<0){
		stop("The number of AR coefficients has to be positive or zero\n")
	}
	return(invisible())

}

validMSM.lm=function(object){
	if(any(object@std<0)) stop("A negative std is not possible!")
	if(length(object@switch)!=ncol(object@Coef)+1){
		stop("The length of sw has to be equal of the number of coefficients in the model plus 1\n")
	}
	return(invisible())

}

validMSM.glm=function(object){
	if(length(object@switch)!=ncol(object@Coef)){
		stop("The length of sw has to be equal of the number of coefficients in the model\n")
	}
	switch(object@family$family,
		poisson=return(invisible()),
		binomial=return(invisible()),
		gaussian=return(invisible()),
		Gamma=return(invisible()),
		stop("The family is not poisson, binomial, gaussian or Gamma!")
	)
	return(invisible())

}




########
#####Class definition
setClass(
	Class= "MSM.linear",
	representation=representation(           
			model  		= "lm",
			Coef  		= "data.frame",
			seCoef 		= "data.frame",
			transMat	= "matrix",
			iniProb 	= "numeric"
	),
	contains="MSM"
)
setClass(
	Class= "MSM.lm",
	representation=representation(           
			std = "numeric"
	),
	contains="MSM.linear"
)

setClass(
	Class= "MSM.glm",
	representation=representation(           
			family = "ANY",
			Likelihood  = "function"
	),
	contains="MSM.linear"
)

####################################
##### Construction, get and set ####
####################################

#### Constructor
msmControl <- function(trace = F,  maxiter = 100, tol = 1e-8, maxiterInner=10, maxiterOuter=5, parallelization=T)
### Control parameters for lmer, glmer and nlmer
{
    stopifnot(maxiter >= 0, tol >= 0, maxiterInner >= 0,maxiterOuter >= 0)
    list(
		maxiter = as.integer(maxiter),
        tol = as.numeric(tol),
		trace = as.logical(trace),
		maxiterInner = as.integer(maxiterInner),
		maxiterOuter = as.integer(maxiterOuter),
		parallelization = as.logical(parallelization)
	)
}


.MSM.formula.msmFit = function(object,k,sw,p,data,family,control){
	call=match.call()
	if (missing(family)) {
		model=lm(object,data) 
	}else{
		model=glm(object,data,family=family) 
	}
	if(missing(p)) p=0 
	data=list(call)
	if(missing(family)&missing(control)){
		msmFit(model,k,sw,p)
	}else{
		if(missing(family)){
			msmFit(model,k,sw,p,control=control)
		}else{
			if(missing(control)){
				msmFit(model,k,sw,p,family=family)
			}else{
				msmFit(model,k,sw,p,family=family,control=control)
			}
		}
	}
}
setMethod(f="msmFit",signature=c("formula","numeric","logical","ANY","data.frame","ANY","ANY"),definition=.MSM.formula.msmFit)


.MSM.lm.msmFit= function(object,k,sw,p,data,family,control){
	if(!missing(data)){
		if(is.list(data)){
			if(class(data[[1]])=="call"){
				call=data[[1]]
			}else{
				call=match.call()
			}
		}else{
			call=match.call()
		}
	}else{
		call=match.call()
	}
	if(missing(p)) p=0 
	if (missing(control)) control=list()
   	control  <- do.call(msmControl, control)

	
	if(p>0){
		var=object$model[,1]
		Ar=apply(as.matrix(1:p),1,function(el){
				length(var)=length(var)-el
				var=c(rep(NA,el),var)
				return(var)	
			}
		)
		colnames(Ar)=paste(names(object$model)[1],"_",1:p,sep="")
		aux=paste(colnames(Ar),collapse="+")
		object=update(formula=as.formula(paste("~.+",aux,sep="")),data=data.frame(object$model,Ar),object)
	
	}

	
	
	Coef=data.frame(matrix(NA,nrow=k,ncol=length(coef(object))))
	std=rep(0,k)

	ind=sample(1:k,length(object$residuals),replace=T)

	for(i in 1:k){
		data1=as.data.frame(object$model[ind==i,,drop=F])
		mod1=update(object,formula=object$terms,data=data1)
		Coef[i,]=coef(mod1)
		std[i]=summary(mod1)$sigma
	}

	names(Coef)=names(coef(object))
	transMat=t(matrix(table(ind,c(ind[-1],NA))/rep(table(ind[-length(ind)]),k),ncol=k))
	ans=new(Class="MSM.lm",
		call=as.call(call),
		model=object,
		k=k,
		switch=sw,
		p=p,
		Coef=Coef,
		std=std,
		transMat=transMat,
		iniProb= rep(1/k,k)
	)
	validMSM.linear(ans)
	validMSM.lm(ans)
	ans=em(ans,control)
	return(ans)
}
setMethod(f="msmFit",signature=c("lm","numeric","logical","ANY","missing","missing","ANY"),definition=.MSM.lm.msmFit)


.MSM.glm.msmFit = function(object,k,sw,p,data,family,control){
	if(!missing(data)){
		if(is.list(data)){
			if(class(data[[1]])=="call"){
				call=data[[1]]
			}else{
				call=match.call()
			}
		}else{
			call=match.call()
		}
	}else{
		call=match.call()
	}
	if(missing(p)) p=0 
	if (missing(control)) control=list()
   	control  <- do.call(msmControl, control)
	
	if(p>0){
		var=object$model[,1]
		Ar=apply(as.matrix(1:p),1,function(el){
				length(var)=length(var)-el
				var=c(rep(NA,el),var)
				return(var)	
			}
		)
		colnames(Ar)=paste(names(object$model)[1],"_",1:p,sep="")
		aux=paste(colnames(Ar),collapse="+")
		object=update(formula=as.formula(paste("~.+",aux,sep="")),data=data.frame(object$model,Ar),object)
	
	}	
	
	
	Coef=data.frame(matrix(NA,nrow=k,ncol=length(coef(object))))
	ind=findInterval(object$residuals,quantile(object$residuals,(1:(k-1))/k))+1
	for(i in 1:k){
		data1=object$model[ind==i,,drop=F]
		mod1=update(object,data=data1)
		Coef[i,]=coef(mod1)
	}
	names(Coef)=names(coef(object))

	transMat=t(matrix(table(ind,c(ind[-1],NA))/rep(table(ind[-length(ind)]),k),ncol=k))


	Likelihood = switch(object$family$family,
			poisson=function(x,mu) dpois(x,lambda=mu),
			binomial=function(x,mu) dbinom(x,prob=mu,size=1),
			gaussian=function(x,mu) dnorm(x,mean=mu,sd=1),
			Gamma=function(x,mu) dgamma(x,shape=mu,rate=1),
			"error"
	)

	ans=new(Class="MSM.glm",
		call=as.call(call),
		model=object,
		k=k,
		switch=sw,
		p=p,
		Coef=Coef,
		transMat=transMat,
		iniProb= rep(1/k,k),

		family=object$family,
		Likelihood=Likelihood
	)
	validMSM.linear(ans)
	validMSM.glm(ans)
	ans=em(ans,control)
	return(ans)
}
setMethod(f="msmFit",signature=c("glm","numeric","logical","ANY","missing","ANY","ANY"),definition=.MSM.glm.msmFit)


##### Get
setMethod(
	f="[",
	signature=c("MSM.lm","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			call = return(x@call),
			model = return(x@model),
			k = return(x@k),
			switch = return(x@switch),
			Coef = return(x@Coef),
			seCoef = return(x@seCoef),
			std = return(x@std),
			transMat = return(x@transMat),
			iniProb = return(x@iniProb),
			Fit = return(x@Fit),
			states = return(x@Fit["states"]),
			stop("Error:",i,"is not a MSM slot")
		)
	}
)


##### Get
setMethod(
	f="[",
	signature=c("MSM.glm","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			call = return(x@call),
			model = return(x@model),
			k = return(x@k),
			switch = return(x@switch),
			Coef = return(x@Coef),
			seCoef = return(x@seCoef),
			transMat = return(x@transMat),
			iniProb = return(x@iniProb),
			Fit = return(x@Fit),
			states = return(x@Fit["states"]),
			family = return(x@family),
			stop("Error:",i,"is not a MSM slot")
		)
	}
)


####################################
########### Functions ##############
####################################

##########
##### show

.MSM.lm.show=function(object){
	cat("Markov Switching Model\n")
	cat("\nCall: ")
	print(object["call"])
	swi=object@switch[-length(object@switch)]
	np=object["k"]*sum(swi)+sum(!swi)
	AIC=2*object["Fit"]["logLikel"]+2*np
	BIC=2*object["Fit"]["logLikel"]+2*np*log(nrow(object@model$model))
	cat("\n")
	print(data.frame(AIC = AIC, BIC = BIC, logLik = -object["Fit"]["logLikel"], 
        row.names = " "))	
	sw=apply(as.matrix(object@switch),1,function(x){
			if(x){
				return("(S)")
			}else{
				return("")
				}
		}
	)
	tau=as.matrix(cbind(object["Coef"],object["std"]))
	dimnames(tau)=list(c(paste(rep("Model",object["k"]),as.character(c(1:object["k"])))),c(paste(c(names(object["Coef"]),"Std"),sw,sep="")))
	cat("\nCoefficients:\n")
   	print(tau)		
	cat("\nTransition probabilities:\n")
	pro=object["transMat"]
	dimnames(pro)=rep(list(paste("Regime",1:object["k"])),2)
	print(pro)
	return(invisible())
}
setMethod(f="show",signature="MSM.lm",definition=.MSM.lm.show)


.MSM.glm.show=function(object){
	cat("Markov Switching Model\n")
	cat("\nCall: ")
	print(object["call"])
	swi=object@switch
	np=object["k"]*sum(swi)+sum(!swi)
	AIC=2*object["Fit"]["logLikel"]+2*np
	BIC=2*object["Fit"]["logLikel"]+2*np*log(nrow(object@model$model))
	cat("\n")
	print(data.frame(AIC = AIC, BIC = BIC, logLik = -object["Fit"]["logLikel"], 
        row.names = " "))	

	sw=apply(as.matrix(object@switch),1,function(x){
			if(x){
				return("(S)")
			}else{
				return("")
			}
		}
	)
	tau=as.matrix(object["Coef"])
	dimnames(tau)=list(c(paste(rep("Model",object["k"]),as.character(c(1:object["k"])))),c(paste(names(object["Coef"]),sw,sep="")))
	print(tau)		
	cat("\nTransition probabilities:\n")
	print(object["transMat"])
	return(invisible())
}
setMethod(f="show",signature="MSM.glm",definition=.MSM.glm.show)



##########
##### summary

.MSM.lm.summary=function(object){
	cat("Markov Switching Model\n")
	cat("\nCall: ")
	print(object["call"])
	swi=object@switch[-length(object@switch)]
	np=object["k"]*sum(swi)+sum(!swi)
	AIC=2*object["Fit"]["logLikel"]+2*np
	BIC=2*object["Fit"]["logLikel"]+2*np*log(nrow(object@model$model))
	cat("\n")
	print(data.frame(AIC = AIC, BIC = BIC, logLik = -object["Fit"]["logLikel"], 
        row.names = " "))	
	sw=apply(as.matrix(swi),1,function(x){
			if(x){
				return("(S)")
			}else{
				return("")
			}
		}
	)
	digits=max(3, getOption("digits") - 3)
		cat("\nCoefficients:\n")
		for(i in 1:object["k"]){
			cat("\nRegime",i,"\n")
			cat("---------\n")
			est=t(round(object["Coef"][i,],digits=digits))
			se=t(round(object["seCoef"][i,],digits=digits))
			nval=round(est/se,digits=digits)
			coefs=cbind(est,se,nval,apply(abs(nval),2,function(el) 2*(1-pnorm(abs(el)))))
			dimnames(coefs) <- list(paste(names(object["Coef"]),sw,sep=""), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
			printCoefmat(coefs)				
			f=object["Fit"]@CondMean[,i]
			r=object["Fit"]@error[,i]
			w=object["Fit"]@smoProb[-1,i]
      		mss <- if (attr(object@model$terms, "intercept")) {
            		m <- sum(w * f/sum(w))
            		sum(w * (f - m)^2)
        			}
		           else sum(w * f^2)
        		rss <- sum(w * r^2)
			r.squared <- mss/(mss + rss)

			cat("\nResidual standard error:",object["std"][i])
			cat("\n")
			cat("Multiple R-squared:", formatC(r.squared, digits = digits))

			resd <- r*sqrt(w)
    			if (length(resd) > 5) {
		        resd <- quantile(resd, na.rm = TRUE)
		        names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
    			}
 			cat("\n\nStandardized Residuals:\n")
			print(resd)		
	}		
	cat("\nTransition probabilities:\n")
	pro=object["transMat"]
	dimnames(pro)=rep(list(paste("Regime",1:object["k"])),2)
	print(pro)
	
	return(invisible())
}
setMethod(f="summary",signature="MSM.lm",definition=.MSM.lm.summary)

.MSM.glm.summary=function(object){
	cat("Markov Switching Model\n")
	cat("\nCall: ")
	print(object["call"])
	swi=object@switch
	np=object["k"]*sum(swi)+sum(!swi)
	AIC=2*object["Fit"]["logLikel"]+2*np
	BIC=2*object["Fit"]["logLikel"]+2*np*log(nrow(object@model$model))
	cat("\n")
	print(data.frame(AIC = AIC, BIC = BIC, logLik = -object["Fit"]["logLikel"], 
        row.names = " "))	
	sw=apply(as.matrix(object@switch),1,function(x){
			if(x){
				return("(S)")
			}else{
				return("")
			}
		}
	)
	digits=max(3, getOption("digits") - 3)
		cat("\nCoefficients:\n")
		for(i in 1:object["k"]){
			cat("\nRegime",i,"\n")
			cat("---------\n")
			est=t(round(object["Coef"][i,],digits=digits))
			se=t(round(object["seCoef"][i,],digits=digits))
			nval=round(est/se,digits=digits)
			coefs=cbind(est,se,nval,apply(abs(nval),2,function(el) 2*(1-pnorm(abs(el)))))
			dimnames(coefs) <- list(paste(names(object["Coef"]),sw,sep=""), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
			printCoefmat(coefs)				


	}		
	cat("\nTransition probabilities:\n")
	pro=object["transMat"]
	dimnames(pro)=rep(list(paste("Regime",1:object["k"])),2)
	print(pro)
	
	return(invisible())
}
setMethod(f="summary",signature="MSM.glm",definition=.MSM.glm.summary)



##########
##### plot

.MSM.plot=function(x,y,...){
	state=apply(x["Fit"]["smoProb"][-1,],1,order,decreasing=T)[1,]
	if((x["k"]+1)>4){
		aux=ifelse(!is.integer((x["k"]+1)/2),1.15+0.15*(x["k"]-4),1)
	}else{
		aux=1
	}
	ts.plot(x["Fit"]["error"][,1],col=2,ylim=c(min(x["Fit"]["error"])-0.5*aux*(abs(max(x["Fit"]["error"]))+abs(min(x["Fit"]["error"])))/2,max(x["Fit"]["error"])),xlab="Order",ylab="Residuals")
	apply(as.matrix(2:x["k"]),1,function(i) 	lines(x["Fit"]["error"][,i],col=1+i))
	lines(apply(as.matrix(1:length(state)),1,function(i) x["Fit"]["error"][i,state[i]]))
	abline(h=0)
	data=x["Fit"]["error"]
	text=apply(as.matrix(1:x["k"]),1,function(i){
			paste("Regime ",i,sep="")
		}
	)
	text[length(text)+1]="Conditional"
	col=apply(as.matrix(1:x["k"]),1,function(i){
			i+1
		}
	)
	col[length(col)+1]=1
	legend(dim(data)[1]*0.5,min(data)-0.1*abs(min(data)),legend=text,cex=0.8,col=col,lty=1,ncol=2)
	return(invisible())
}
setMethod(f="plot",signature=c("MSM.linear","missing"),definition=.MSM.plot)


##########
##### plotProb

.MSM.plotProb=function(x,which){
	if(missing(which)){
		which=1:(x["k"]+1)
	}else{
		if(is.numeric(which)){
			if(any(which > (x["k"]+1))|any(which < 1)) stop("You must to write numbers between 1 and 3.")
		}else{
			stop("You must write numbers.")
		}
	}
  oldpar=par()
	a=ceiling(x["k"]/3)
	aux1=1:a
	aux2=matrix(NA,nrow=x["k"],ncol=1)
	aux2[,1]=apply(as.matrix(1:x["k"]),1,function(i){
			if(any((aux1*3)==(i-1))){
				T
			}else{
				F
			}
		}
	)
	aux2[1,1]=T
	cont=1
	if(any(which==1)){	
		for (i in 1:x["k"]){		
			if(aux2[i,1]){
				if(i>1){
					if(cont<length(aux1)){
						par(omi=c(0.1,0.1,0.1,0.5))
						par(mfrow=c(3,1))
						cont=cont+1
					}else{
						par(omi=c(0.1,0.1,0.1,0.5))
						par(mfrow=c(x["k"]-(a-1)*3,1))
					}
				}else{
					if(x["k"]/3<1){
						par(omi=c(0.1,0.1,0.1,0.5))
						par(mfrow=c(2,1))
					}else{
						par(omi=c(0.1,0.1,0.1,0.5))
						par(mfrow=c(3,1))
					}
					cont=cont+1
				}
			}
			plot(x["Fit"]["filtProb"][,i],main=paste("Regime",i),ylim=c(0,1),xlab="",ylab="",type="h")
			lines(x["Fit"]["smoProb"][-1,i],col=2)
			par(las=3)
			mtext("Smoothed Probabilities",side=2,line=2.5,col=2)
			mtext("Filtered Probabilities",side=4,line=2.5)
		}
	}
	
	if(any(which>1)){	
		aux=which[which>1]-1
		z=x["model"]$model[1]	
		apply(as.matrix(1:length(aux)),1,function(i){
				a=layout(matrix(c(1,1,1,2),ncol=1,nrow=4),TRUE)
				y=x["Fit"]["smoProb"][-1,aux[i]]
				par(omi=c(0.1,0.1,0.1,0.5))
				par(las=1,yaxt="n")
				plot(0,type="l",xlim=c(1,length(t(z))),ylim=c(min(z),max(z)),main=paste("Regime",aux[i]),xlab=paste(names(z),"vs. Smooth Probabilities"),ylab="")
				val=cbind(which(diff(c(0,findInterval(y,0.5)))==1),which(diff(c(findInterval(y,0.5),0))==-1))
				apply(val,1,function(el) rect(el[1],min(z),el[2],max(z),col="light grey",border=NA))
				par(new=T,las=1,bty="o",yaxt="n")			
				plot(ts(z),col=1,ylim=c(min(z),max(z)),xlab="",ylab="")
				par(las=3,yaxt="s")
				mtext(names(z),side=2,line=2.5,col=1)
				axis(side=4)
				barplot(x["Fit"]["smoProb"][-1,aux[i]],ylim=c(0,1))	
			}
		)
	}
  par(mfrow=c(1,1))
	return(invisible())
}
setMethod(f="plotProb",signature=c("MSM.linear","ANY"),definition=.MSM.plotProb)



##########
##### plotReg

.MSM.larg.plotReg=function(x,expl,regime){
	if(missing(regime)){
		regime=1
	}else{
		if(is.numeric(regime)){
			if(any(regime > x["k"])|any(regime < 1)) stop("You must to write a correct regime.")
		}else{
			stop("regime must to be numeric.")
		}
	}
	z=x["model"]$model[1]
	apply(as.matrix(1:length(regime)),1,function(j){
		apply(as.matrix(1:length(x["model"]$model[-1])),1,function(i){
				a=layout(matrix(c(1,2,1,2),ncol=1,nrow=2),TRUE)
				y=x["Fit"]["smoProb"][-1,regime[j]]
				v=x["model"]$model[i+1]
				par(omi=c(0.1,0.01,0.1,0.1))
				par(las=1)
				plot(0,type="l",xlim=c(1,length(t(z))),ylim=c(min(z),max(z)),main=paste("Regime",regime[j],sep=""),xlab=paste(names(z),"with Smooth Probabilities"),ylab="")
				val=cbind(which(diff(c(0,findInterval(y,0.5)))==1),which(diff(c(findInterval(y,0.5),0))==-1))
				apply(val,1,function(el) rect(el[1],min(z),el[2],max(z),col="light grey",border=NA))
				par(new=T,las=1,bty="o")
				plot(ts(z),col=1,ylim=c(min(z),max(z)),main="",xlab="",ylab="")
				plot(0,type="l",xlim=c(1,length(t(v))),ylim=c(min(v),max(v)),main="",xlab=paste(names(v),"with Smooth Probabilities"),ylab="")
				apply(val,1,function(el) rect(el[1],min(v),el[2],max(v),col="light grey",border=NA))
				par(new=T,las=1,bty="o")
				plot(ts(v),col=1,ylim=c(min(v),max(v)),main="",xlab="",ylab="")
			}
			)
		}
	)
	par(mfrow=c(1,1))
	return(invisible())
}
setMethod(f="plotReg",signature=c("MSM.linear","missing","ANY"),definition=.MSM.larg.plotReg)


.MSM.sma.plotReg=function(x,expl,regime){
	apply(as.matrix(expl),1,function(var) if(!any(apply(as.matrix(names(x["model"]$model[-1])),1,function(el) ifelse(var==el,T,F)))) stop("The name of the variable is not correct.\n"))
	if(missing(regime)){
		regime=1
	}else{
		if(is.numeric(regime)){
			if(any(regime > x["k"])|any(regime < 1)) stop("You must to write a correct regime.")
		}else{
			stop("regime must to be numeric.")
		}
	}
	z=x["model"]$model[1]
	apply(as.matrix(1:length(regime)),1,function(j){
		apply(as.matrix(1:length(expl)),1,function(i){
				a=layout(matrix(c(1,2,1,2),ncol=1,nrow=2),TRUE)
				y=x["Fit"]["smoProb"][-1,regime[j]]
				v=x["model"]$model[expl[i]]
				par(omi=c(0.1,0.1,0.1,0.1))
				par(las=1,yaxt="n")
				plot(0,type="l",xlim=c(1,length(t(z))),ylim=c(min(z),max(z)),main=paste("Regime",regime[j],sep=""),xlab=paste(names(z),"with Smooth Probabilities"),ylab="")
				val=cbind(which(diff(c(0,findInterval(y,0.5)))==1),which(diff(c(findInterval(y,0.5),0))==-1))
				apply(val,1,function(el) rect(el[1],min(z),el[2],max(z),col="light grey",border=NA))
				par(new=T,las=1,bty="o",yaxt="n")
				plot(ts(z),col=1,ylim=c(min(z),max(z)),main="",xlab="",ylab="")
				plot(0,type="l",xlim=c(1,length(t(v))),ylim=c(min(v),max(v)),main="",xlab=paste(names(v),"with Smooth Probabilities"),ylab="")
				apply(val,1,function(el) rect(el[1],min(v),el[2],max(v),col="light grey",border=NA))
				par(new=T,las=1,bty="o",yaxt="n")
				plot(ts(v),col=1,ylim=c(min(v),max(v)),main="",xlab="",ylab="")
			}
			)
		}
	)
	par(mfrow=c(1,1))
	return(invisible())
}
setMethod(f="plotReg",signature=c("MSM.linear","character","ANY"),definition=.MSM.sma.plotReg)

##########
##### plotDiag

.MSM.larg.plotDiag=function(x,regime,which){
	if(missing(which)){
		which=1:3
	}else{
		if(is.numeric(which)){
			if(any(which > 3)|any(which < 1)) stop("You must to write numbers between 1 and 3.")
		}else{
			stop("You must write numbers.")
		}
	}
	residPooled=apply(x["Fit"]["error"]*x["Fit"]["smoProb"][-1,],1,sum)
	if(any(which==1)){	
		ts.plot(residPooled,main="Pooled residuals")
		abline(h=0)
		abline(h=c(-3*sd(residPooled),3*sd(residPooled)),lty=3,col=4)
	}
	if(any(which==2)){
		qqnorm(residPooled,main="Normal Q-Q Plot Pooled Residuals")
		qqline(residPooled,col=2,lwd=2)
	}
	if(any(which==3)){
		par(mfrow=c(2,2))
		acf(residPooled,ylim=c(-1,1),main="ACF of Residuals")
		pacf(residPooled,ylim=c(-1,1),main="PACF of Residuals")
		acf(residPooled^2,ylim=c(-1,1),main="ACF of Square Residuals")
		pacf(residPooled^2,ylim=c(-1,1),main="PACF of Square Residuals")
		par(mfrow=c(1,1))
	}
	return(invisible())
}
setMethod(f="plotDiag",signature=c("MSM.linear","missing","ANY"),definition=.MSM.larg.plotDiag)

.MSM.larg.plotDiag=function(x,regime,which){
	if(regime[length(regime)]=="all"){
		aux=c(1:x["k"])
	}else{
		if(any(regime > x["k"])|any(regime < 1)) stop("You must to write a correct regime.")
		aux=regime
	}
	if(missing(which)){
		which=1:3
	}else{
		if(is.numeric(which)){
			if(any(which > 3)|any(which < 1)) stop("You must to write numbers between 1 and 3.")
		}else{
			stop("You must write numbers.")
		}
	}
	apply(as.matrix(aux),1,function(i){
			if(any(which==1)){
				ts.plot(x["Fit"]["error"][,i],main=paste("Regime ",i,sep=""),ylab="Residuals")
				abline(h=0)
				abline(h=c(-3*sd(x["Fit"]["error"][,i]),3*sd(x["Fit"]["error"][,i])),lty=3,col=4)
			}
			if(any(which==2)){
				qqnorm(x["Fit"]["error"][,i],main=paste("Normal Q-Q Plot Regime ",i,sep=""))
				qqline(x["Fit"]["error"][,i],col=2,lwd=2)
			}
			if(any(which==3)){			
				par(mfrow=c(2,2))
				acf(x["Fit"]["error"][,i],ylim=c(-1,1),main=paste("ACF of Residuals. Reg: ",i,sep=""))
				pacf(x["Fit"]["error"][,i],ylim=c(-1,1),main=paste("PACF of Residuals. Reg: ",i,sep=""))
				acf(x["Fit"]["error"][,i]^2,ylim=c(-1,1),main=paste("ACF of Square Resid. Reg: ",i,sep=""))
				pacf(x["Fit"]["error"][,i]^2,ylim=c(-1,1),main=paste("PACF of Square Resid. Reg: ",i,sep=""))
				par(mfrow=c(1,1))
			}
		}
	)
	return(invisible())
}
setMethod(f="plotDiag",signature=c("MSM.linear","ANY","ANY"),definition=.MSM.larg.plotDiag)



##########
##### resid
.MSM.lm.sma.msmResid=function(object,regime){
	return(apply(object["Fit"]["error"]*object["Fit"]["smoProb"][-1,],1,sum))
}
setMethod(f="msmResid",signature=c("MSM.lm","missing"),definition=.MSM.lm.sma.msmResid)

.MSM.lm.larg.msmResid=function(object,regime){
	if(regime[length(regime)]=="all"){
		aux=c(1:object["k"])
		res=apply(as.matrix(aux),1,function(i) object["Fit"]["error"][,i])
		dimnames(res)[[2]]=paste("Regime ",1:object["k"],sep="")
		return(res)		
	}else{
		if(any(regime > object["k"])|any(regime < 1)) stop("You must to write a correct regime.")
		aux=regime
		return(object["Fit"]["error"][,regime])
		}
}
setMethod(f="msmResid",signature=c("MSM.lm","ANY"),definition=.MSM.lm.larg.msmResid)

.MSM.glm.sma.msmResid=function(object,regime){
	res=object["Fit"]["error"]/sqrt(object@family$variance(object["Fit"]@CondMean))
	return(apply(res*object["Fit"]["smoProb"][-1,],1,sum))
}
setMethod(f="msmResid",signature=c("MSM.glm","missing"),definition=.MSM.glm.sma.msmResid)

.MSM.glm.larg.msmResid=function(object,regime){
	if(regime[length(regime)]=="all"){
		aux=c(1:object["k"])
		res=apply(as.matrix(aux),1,function(i) object["Fit"]["error"][,i]/sqrt(object@family$variance(object["Fit"]@CondMean[,i])))
		dimnames(res)[[2]]=paste("Regime ",1:object["k"],sep="")
		return(res)		
	}else{
		if(any(regime > object["k"])|any(regime < 1)) stop("You must to write a correct regime.")
		aux=regime
		return(object["Fit"]["error"][,regime]/sqrt(object@family$variance(object["Fit"]@CondMean[,regime])))
		}
}
setMethod(f="msmResid",signature=c("MSM.glm","ANY"),definition=.MSM.glm.larg.msmResid)


AIC.MSM.lm <-
  
  function(object, ..., k=2)
{
	swi=object@switch
	np=object["k"]*sum(swi)+sum(!swi)
	return(2*object["Fit"]["logLikel"]+k*np)
}
AIC.MSM.glm <-
  
  function(object, ..., k=2)
{
	swi=object@switch
	np=object["k"]*sum(swi)+sum(!swi)
	return(2*object["Fit"]["logLikel"]+k*np)
}

AIC <-
  ## Return the object's value of the Bayesian Information Criterion
  function(object, ...,k=2) UseMethod("AIC")

intervals.MSM.lm=function(object,level=0.95,...){
	cat("\nAproximate intervals for the coefficients. Level=",level,"\n")
	aux=names(object["Coef"])
	lower=object["Coef"]-qnorm(1-(1-level)/2)*object["seCoef"]
	upper=object["Coef"]+qnorm(1-(1-level)/2)*object["seCoef"]
	a=apply(as.matrix(1:length(aux)),1,function(i){
			cat(paste("\n",aux[i],": \n",sep=""))
			#cat("---------\n")
			intmat=cbind(lower[aux[i]],object["Coef"][aux[i]],upper[aux[i]])
			dimnames(intmat)=list(c(paste("Regime ",1:object["k"],sep="")),c("Lower","Estimation","Upper"))
			print(intmat)
			#cat("---------\n")
			cat("\n")
		}
	)
}
intervals.MSM.glm=function(object,level=0.95,...){
	cat("\nAproximate intervals for the coefficients. Level=",level,"\n")
	aux=names(object["Coef"])
	lower=object["Coef"]-qnorm(1-(1-level)/2)*object["seCoef"]
	upper=object["Coef"]+qnorm(1-(1-level)/2)*object["seCoef"]
	a=apply(as.matrix(1:length(aux)),1,function(i){
			cat(paste("\n",aux[i],": \n",sep=""))
			#cat("---------\n")
			intmat=cbind(lower[aux[i]],object["Coef"][aux[i]],upper[aux[i]])
			dimnames(intmat)=list(c(paste("Regime ",1:object["k"],sep="")),c("Lower","Estimation","Upper"))
			print(intmat)
			#cat("---------\n")
			cat("\n")
		}
	)
}
intervals <-
  ## Return the object's value of the Bayesian Information Criterion
  function(object,level=0.95,...) UseMethod("intervals")


#########
#### msmfilter

.MSM.lm.msmFilter= function(object){
	model=object["model"]
	k=object["k"]
	
	Coef=object["Coef"]
	std=object["std"]
	P=object["transMat"]

	# Calculation of some preliminar variables
	nr=length(model$model[,1])
	terms=model.matrix(model)

	CondMean=as.matrix(terms)%*%t(as.matrix(Coef))
	error= as.matrix(model$model[,1,drop=F])%*%matrix(rep(1,k),nrow=1)-CondMean
	Likel=t(dnorm(t(error),0,std))

	###Filtered Probabilities ####
	fProb=matrix(data=0,nrow=nr,ncol=k)
	margLik=matrix(data=0,nrow=nr,ncol=1)
		
	fProb[1,]= (P%*%matrix(object["iniProb"],ncol=1))*t(Likel[1,,drop=F])
	margLik[1,1] = sum(fProb[1,])
	fProb[1,] = fProb[1,] / margLik[1,1]
	for (i in 2:nr){
		# Mixtura de funcions
		# MS filter equation
		# MS filter Filter margLikuation for probabilities
		fProb[i,] = (P%*%t(fProb[i-1,,drop=F]))*t(Likel[i,,drop=F])
		margLik[i,1] = sum(fProb[i,])
		fProb[i,] = fProb[i,]/margLik[i,1]
	}

	# Negative sum of log Likelihood 
	loglik=-sum(log(margLik[1:nr]))

	# Passing up to output structure
	ans=new(Class="MSM.fitted",CondMean=CondMean,error=error, Likel=Likel,margLik=margLik, filtProb=fProb, logLikel=loglik )
	return(ans)
}
setMethod(f="msmFilter",signature=c("MSM.lm"),definition=.MSM.lm.msmFilter)

.MSM.glm.msmFilter=function(object){
	model=object["model"]
	k=object["k"]
	family=model$family

	Coef=object["Coef"]
	P=object["transMat"]

	# Calculation of some preliminar variables
	nr=length(model$model[,1])

	terms=model.matrix(model)
	
	CondMean=family$linkinv(as.matrix(terms)%*%t(as.matrix(Coef)))
	error= as.matrix(model$model[,1,drop=F])%*%matrix(rep(1,k),nrow=1)-CondMean

	Likel=object@Likelihood(as.matrix(model$model[,1,drop=F])%*%matrix(rep(1,k),nrow=1),mu=CondMean)

	###Filtered Probabilities ####
	fProb=matrix(data=0,nrow=nr,ncol=k)
	margLik=matrix(data=0,nrow=nr,ncol=1)
	
	margLik[1,1]=sum ((P%*%matrix(object["iniProb"],ncol=1)) * t(Likel[1,,drop=F]))
	fProb[1,]= ((P%*%matrix(object["iniProb"],ncol=1))*t(Likel[1,,drop=F]))/margLik[1,1]
	for (i in 2:nr){
		# Mixtura de funcions
		# MS filter margLikuation
		margLik[i,1]=sum ((P%*%t(fProb[i-1,,drop=F])) * t(Likel[i,,drop=F]))
		# MS filter Filter margLikuation for probabilities
		fProb[i,]= ((P%*%t(fProb[i-1,,drop=F])*t(Likel[i,,drop=F]))/margLik[i,1])
	}

	# Negative sum of log Likelihood for fmincon (fmincon minimzes the function)
	loglik=-sum(log(margLik[1:nr]))

	# Passing up to output structure
	ans=new(Class="MSM.fitted",CondMean=CondMean,error=error, Likel=Likel,margLik=margLik, filtProb=fProb, logLikel=loglik )
	return(ans)
}
setMethod(f="msmFilter",signature=c("MSM.glm"),definition=.MSM.glm.msmFilter)

#########
#### msmsmooth

.MSM.msmSmooth=function(object){		
	object@Fit=msmFilter(object)
	nr=length(object["model"]$model[,1])
	fProb=object["Fit"]["filtProb"]
	k=object["k"]
	P=object["transMat"]
	smoTransMatrob=matrix(0,ncol=k,nrow=nr+1)
	smoTransMatrob[nr+1,]=fProb[nr,]
	#smoTransMatrob[nr+1,]=fProb[nr,]%*%t(P)

	proba=rbind(object@iniProb,fProb)
	pro=proba%*%t(P)
	smoTransMat=list(NULL)
	for (i in (nr-1):0){
		smoTransMat[[i+1]]=matrix(0,ncol=k,nrow=k)
		for (ini in 1:k){
			smoTransMatrob[i+1,ini]=0
			for (fi in 1:k){
				smoTransMat[[i+1]][ini,fi]=smoTransMatrob[i+2,fi]*proba[i+1,ini]*P[fi,ini]/pro[i+1,fi]
				smoTransMatrob[i+1,ini]=smoTransMatrob[i+1,ini] + smoTransMat[[i+1]][ini,fi]
			}
		}
	}
	object@Fit@smoProb=smoTransMatrob
	object@Fit@smoTransMat=smoTransMat
	return(object)
}
setMethod(f="msmSmooth",signature=c("MSM.linear"),definition=.MSM.msmSmooth)


#########
####optimizer
fopt.lm=function(param, object=object){
	if (tail(object["switch"],1)==F){
		object@std<-exp(param[1])
		ini=1
	} else {
		object@std<-exp(param[1:object["k"]])
		ini=object["k"]
 	}
	long=ini+(object["k"]-1)*object["k"]
	mprob=matrix(logitinv(c(param[(ini+1):long])),nrow=object["k"],byrow=T)
	object@transMat<-matrix(c(mprob,1-apply(mprob,1, function(x) sum(x))),nrow=object["k"],byrow=T)
	swi=object["switch"][-length(object["switch"])]
	mi=sum(!swi)
	aux=object["Coef"]
	aux[,which(swi)]=as.data.frame(matrix(param[-c(1:(long+mi))],nrow=object["k"],byrow=T))
	aux[,which(!swi)]=as.data.frame(matrix(rep(param[long+(1:mi)],object["k"]),nrow=object["k"],byrow=T))
	object@Coef=aux
	return(msmFilter(object)@logLikel)
}
fopt.glm=function(param, object=object){
	long=(object["k"]-1)*object["k"]
	mprob=matrix(logitinv(c(param[1:long])),ncol=object["k"],byrow=T)
	object@transMat<-matrix(c(mprob,1-apply(mprob,2, function(x) sum(x))),nrow=object["k"],byrow=T)
	swi=object["switch"]	
	mi=sum(!swi)
	aux=object["Coef"]
	aux[,which(swi)]=as.data.frame(matrix(param[-c(1:(long+mi))],nrow=object["k"],byrow=T))
	aux[,which(!swi)]=as.data.frame(matrix(rep(param[long+(1:mi)],object["k"]),nrow=object["k"],byrow=T))
	object@Coef=aux
	return(msmFilter(object)@logLikel)
}


.MSM.lm.hessian=function(object){
		if (tail(object["switch"],1)==F){
			lstd=log(object["std"][1])
		} else {
			lstd=log(object["std"])
	 	} 
		swi=object["switch"][-length(object["switch"])]
   		param=c(lstd,
		logit(matrix(object["transMat"][1:object["k"]-1,],nrow=1,byrow=T)),
		object["Coef"][1,!swi], 
		matrix(t(as.matrix(object["Coef"])[,swi]),nrow=1))
		res=fdHess(
			pars=param,
			fun=fopt.lm,
			object=object
		)
		long=object["k"]+(object["k"]-1)*object["k"]
		mi=sum(!swi)
	      hessian=sqrt(abs(diag(solve(res$Hessian))))
		stdaux=object["Coef"]
		stdaux[,which(swi)]=as.data.frame(matrix(hessian[-c(1:(long+mi))],nrow=object["k"],byrow=T))
		stdaux[,which(!swi)]=as.data.frame(matrix(rep(hessian[long+(1:mi)],object["k"]),nrow=object["k"],byrow=T))
		object@seCoef=stdaux
		return(object)  
}
setMethod(f="hessian",signature=c("MSM.lm"),definition=.MSM.lm.hessian)


.MSM.glm.hessian=function(object){
    		param=c(
		logit(matrix(object["transMat"][1:object["k"]-1,],nrow=1,byrow=T)),
		as.matrix(object["Coef"])[1,!object["switch"]], 
		matrix(t(as.matrix(object["Coef"])[,object["switch"]]),nrow=1))
		res=fdHess(
			pars=param,
			fun=fopt.glm,
			object=object
		)
		long=(object["k"]-1)*object["k"]
		mi=sum(!object["switch"])
	      hessian=sqrt(abs(diag(solve(res$Hessian))))
		stdaux=object["Coef"]
		stdaux[,which(object["switch"])]=as.data.frame(matrix(hessian[-c(1:(long+mi))],nrow=object["k"],byrow=T))
		stdaux[,which(!object["switch"])]=as.data.frame(matrix(rep(hessian[long+(1:mi)],object["k"]),nrow=object["k"],byrow=T))
		object@seCoef=stdaux
		return(object)  
}
setMethod(f="hessian",signature=c("MSM.glm"),definition=.MSM.glm.hessian)

.MSM.lm.maximEM=function(object,dades){
 	k=object["k"]
	swi=object["switch"][-length(object["switch"])]
	co=object["Coef"]
	w=object["Fit"]["smoProb"][-1,]

	modaux=lm(y~.-1,dades,weights=c(t(w)))	
	object@Coef=as.data.frame(matrix(rep(coef(modaux),rep(ifelse(swi,1,k),ifelse(swi,k,1))),nrow=k))
	if (tail(object["switch"],1)==T){
		object@std=sqrt(apply(w*matrix(resid(modaux),ncol=k,byrow=T)^2,2,sum)/apply(w,2,sum))
	} else {
		std=sum(weighted.residuals(modaux)^2)/nrow(w)
		object@std=rep(sqrt(std),k)
	}
	names(object@Coef)=names(co)
	return(object)
}

setMethod(f="maximEM",signature=c("MSM.lm","data.frame"),definition=.MSM.lm.maximEM)

.MSM.glm.maximEM=function(object,dades){
 	k=object["k"]
	sw=object["switch"]
	co=object["Coef"]
	w=object["Fit"]["smoProb"][-1,]

	modaux=glm(y~.-1,dades,weights=c(t(w)),family=object["model"]$family)	
	
	object@Coef=as.data.frame(matrix(rep(coef(modaux),rep(ifelse(sw,1,k),ifelse(sw,k,1))),nrow=k))
	names(object@Coef)=names(co)
	return(object)

}
setMethod(f="maximEM",signature=c("MSM.glm","data.frame"),definition=.MSM.glm.maximEM)


.MSM.linear.iteraEM=function(object,dades,control){
	k=object["k"]
	for (it in 1:control$maxiter){
		oldcoef=object["Coef"]
		oldll=object["Fit"]["logLikel"]

		##M-step
		object=maximEM(object,dades)
		smoTransMatrob=object["Fit"]["smoProb"]
		smoTransMat=object["Fit"]["smoTransMat"]
		object@transMat=matrix(apply(matrix(unlist(smoTransMat),nrow=k*k),1,sum)/rep(apply(smoTransMatrob[-1,],2,sum),rep(k,k)),ncol=k)
		object@iniProb=object["Fit"]["smoProb"][1,]

		##E-step
		object=msmSmooth(object)

		if (control$trace) cat(" Inner Iter.",it," logLikel=",object["Fit"]["logLikel"],"\n")
		if ((max(abs(object["Fit"]["logLikel"] - oldll))/(0.1 + max(abs(object["Fit"]["logLikel"]))) < control$tol)& (max(abs(object["Coef"] - oldcoef))/(0.1 + max(abs(object["Coef"]))) < control$tol)) break
	}
	return(object)
}
setMethod(f="iteraEM",signature=c("MSM.linear","data.frame","ANY"),definition=.MSM.linear.iteraEM)



####em
.MSM.em=function(object,control){
	
  		k=object["k"]
		swi=object["switch"]
		co=object["Coef"]
		
		constX=function(el,swit){
			if(swit){ mat=diag(1,k) }else{ mat=rep(1,k)}
			kronecker(el,mat)
		}
		Xini=model.matrix(object["model"])
		X=NULL
		for(i in 1:ncol(Xini)){ X=cbind(X,constX(Xini[,i,drop=F],swi[i])) }
		y=kronecker(as.matrix(object["model"]$model[,1,drop=F]),rep(1,k))
		dades=data.frame(y=y,X)		
		
		mprob=object@transMat[-1,,drop=F]
		object@transMat=matrix(c(1-apply(mprob,2, function(x) sum(x)),mprob),nrow=object["k"],byrow=T)
		object@transMat[object@transMat<0]=0
		object@transMat[object@transMat>1]=1
		object=msmSmooth(object)
		
		maxiterInner=control$maxiterInner
		maxiterOuter=control$maxiterOuter
		parallelization=control$parallelization
      
  	if(parallelization){
  		    mc=detectCores(logical = TRUE)
  		    cl <- makeCluster(mc)
  		    #clusterExport(cl,c("dades","object","control","maxiterInner"))
  	}
	
		paralel=function(id){
			x<-object
			smoTransMat=lapply(vector("list",nrow(x@Fit@filtProb)),function(el){
					ma=runif(k*k)
			   	 	matrix(ma/sum(ma),ncol=k)
				})
			smoTransMatrob=rbind(t(sapply(smoTransMat,function(el)apply(el,1,sum))),runif(k))
			
			x@Fit@smoTransMat<-smoTransMat
			x@Fit@smoProb=smoTransMatrob
			x=iteraEM(x,dades,control=list(maxiter=maxiterInner,tol=control$tol,trace=control$trace,parallelization=control$parallelization))
			return(list(Minim=x@Fit@logLikel,inismoTransMat=x@Fit@smoTransMat,inismoTransMatrob=x@Fit@smoProb))
		}
  	if(parallelization){
  	      junk <- clusterEvalQ(cl, library(MSwM))
					paralRes=parLapply(cl, c(1:maxiterOuter),paralel)
          stopCluster(cl)
  		  } else {
  		    paralRes=lapply(c(1:maxiterOuter),paralel)
  		  }
        
				
		
		Minim=paralRes[[1]][["Minim"]]
		inismoTransMat=paralRes[[1]][["inismoTransMat"]]
		inismoTransMatrob=paralRes[[1]][["inismoTransMatrob"]]
		i=2
		while(i<=maxiterOuter){
			if(Minim>paralRes[[i]][["Minim"]]){
				Minim=paralRes[[i]][["Minim"]]
				inismoTransMat=paralRes[[i]][["inismoTransMat"]]
				inismoTransMatrob=paralRes[[i]][["inismoTransMatrob"]]
			}
			i=i+1
		}

		if (control$trace) cat("Initial Value:",Minim,"\n")	
  		object@Fit@smoTransMat=inismoTransMat
  		object@Fit@smoProb=inismoTransMatrob
		object=iteraEM(object,dades,control)
		if (control$trace) cat("Calculating standard errors...\n")
		object=hessian(object)
		return(object)
	
}
setMethod(f="em",signature=c("MSM.linear","list"),definition=.MSM.em)


