##########################################################################
## Show and plot methods
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

##*****************************************************************************
##                        Print Method
##*****************************************************************************          

show.SAVE<- function(object){
	result<- list()	
#only STAGEI:
	result$stage2<- FALSE
	result$call<- object@call
	result$mle<- object@mle
#STAGEII also:
	if (dim(object@mcmcsample)[1]>0){
		result$bayesfitcall<- object@bayesfitcall
		result$stage2<- TRUE
		mysummary<- function(x){
			round(c(mean(x),sd(x),median(x),
					quantile(x,probs=c(.025,.25,.75,.975))),3)
		}
		result$mcmcsum<- matrix(0, nrow=dim(object@mcmcsample)[2], ncol=7)
		result$mcmcsum<- t(apply(object@mcmcsample, MARGIN=2, FUN=mysummary))
		rownames(result$mcmcsum)<- colnames(object@mcmcsample)
		colnames(result$mcmcsum)<- c("Mean", "SD", "Median", "2.5%", 
									 "25%", "75%", "97.5%")
	}
	cat("\n")
	cat("---------------\n")
	cat("call to SAVE:\n")
	print(result$call)
	cat("---------------\n")
	cat("Maximum likelihood estimates:\n")
	print(result$mle)
	if(!result$stage2){cat("---------------\n Message:bayesfit has not been run.\n")}
	if(result$stage2){
		cat("---------------\n")
		cat("call to bayesfit:\n")		
		print(result$bayesfitcall)
		cat("---------------\n")
		cat("Summaries of the MCMC:\n")
		print(result$mcmcsum)
	}
	cat("\n\n\n")
}

if(!isGeneric("show")) {
        #setGeneric(name = "show", simpleInheritanceOnly = TRUE)
	setGeneric(name = "show",
             def = function(object) standardGeneric("show")
	)
}

setMethod ("show", "SAVE", function(object){
	show.SAVE (object)
	}
)


##*****************************************************************************
##                        Plot Method
##*****************************************************************************
plot.SAVE<- function(x, option="trace",...){
		object <- x
		if (dim(object@mcmcsample)[1]==0){stop("Nothing to be plotted before bayesfit is run.\n")}

		if (option=="trace")
		{
			target<- as.mcmc(object@mcmcsample)
			plot(target, trace=TRUE, density=FALSE, ...)
		}
		if (option=="calibration"){
			if (length(object@calibrationnames)!=0)
			{
			 howmanycal<- dim(object@mcmcsample)[2]-2
			 nc<- round(sqrt(howmanycal))
			 nr<- ceiling(howmanycal/nc)
			 par(mfrow=c(nr,nc))
			 for (i in 1:howmanycal){

				#thisprior<- as.numeric(object@prior[(i-1)*6+2:6])
				#thisprior<- as.numeric(object@prior[(i-1)*5+1:5])
				thisprior<- as.numeric(object@prior[i,1:5])

				hs<- hist(object@mcmcsample[,i], breaks=seq(from=thisprior[2], to=thisprior[3], length=20), 
						  col=gray(.85),border=gray(1), main="", 
						  xlab=colnames(object@mcmcsample)[i], prob=T)

				xpr<- seq(from=thisprior[2],to=thisprior[3],length=100)
				if (thisprior[1]==1){
					densprior<- function(x){
						denom<- pnorm(q=thisprior[3], mean=thisprior[4], sd=sqrt(thisprior[5]))-pnorm(q=thisprior[2], mean=thisprior[4], sd=sqrt(thisprior[5]))
						dnorm(x, mean=thisprior[4], sd=sqrt(thisprior[5]))/denom
					}}
				if (thisprior[1]==0){
					densprior<- function(x){
						dunif(x, min=thisprior[2], max=thisprior[3])	
					}}

				ypr<- vapply(xpr, FUN=densprior, FUN.VALUE=1)
				lines(xpr, ypr, type="l", ...)				

			 }
			} else print("There are no calibration parameters. Nothing to plot.")
		}
		if (option=="precision")
		{
			par(mfrow=c(1,2))
			ncal<- dim(object@mcmcsample)[2]-2
			#plot lambdaB
			hs<- hist(object@mcmcsample[,ncal+1], col=gray(.85),border=gray(1), main="", 
				 xlab=colnames(object@mcmcsample)[ncal+1], prob=T, breaks=20)
			xpr<- seq(from=hs$breaks[1], to=hs$breaks[length(hs$breaks)],length=100)
			ypr<- dexp(x=xpr,rate=1/(object@mle$thetaF['lambdaB']*object@mcmcMultmle))
			lines(xpr, ypr, type="l", ...)
			abline(v=object@mle$thetaF['lambdaB'],lty=2)
			#print(object@mle$thetaF['lambdaB'])
			
			#plot lambdaF
			hs<- hist(object@mcmcsample[,ncal+2], col=gray(.85),border=gray(1), main="", 
				 xlab=colnames(object@mcmcsample)[ncal+2], prob=T, breaks=20)
			xpr<- seq(from=hs$breaks[1], to=hs$breaks[length(hs$breaks)],length=100)
			ypr<- dexp(x=xpr,rate=1/(object@mle$thetaF['lambdaF']*object@mcmcMultmle))
			lines(xpr, ypr, type="l", ...)
			#print(object@mle$thetaF['lambdaF'])
			abline(v=object@mle$thetaF['lambdaF'],lty=2)
		}
	par(mfrow=c(1,1))
	
}

#if(!isGeneric("plot")) {
#	setGeneric(name = "plot",
#	def = function(x, y, ...) standardGeneric("plot")
#	)
#}

setMethod(f="plot",
	signature(x = "SAVE",y="missing"), 
	function(x, option="trace", ...) {
	plot.SAVE(x = x, option = option, ...)
}
)

##*****************************************************************************
##                        Summary
##*****************************************************************************

summary.SAVE<- function(object,...){
	result<- new("summary.SAVE")	
#only STAGEI:
	result@stage2<- FALSE
	result@call<- object@call
	result@mle<- object@mle
#STAGEII also:
	if (dim(object@mcmcsample)[1]>0){
		result@bayesfitcall<- object@bayesfitcall
		result@stage2<- TRUE
		mysummary<- function(x){
			round(c(mean(x),sd(x),median(x),
					quantile(x,probs=c(.025,.25,.75,.975))),3)
		}
		result@mcmcsum<- matrix(0, nrow=dim(object@mcmcsample)[2], ncol=7)
		result@mcmcsum<- t(apply(object@mcmcsample, MARGIN=2, FUN=mysummary))
		rownames(result@mcmcsum)<- colnames(object@mcmcsample)
		colnames(result@mcmcsum)<- c("Mean", "SD", "Median", "2.5%", 
									 "25%", "75%", "97.5%")
	}
	return(result)
}


if(!isGeneric("summary")) {
	setGeneric(name = "summary",
			   def = function(object, ...) standardGeneric("summary")
			   )
}

setMethod("summary","SAVE",
	function(object){ summary.SAVE (object) }
)



show.summary.SAVE<- function(object){
	cat("\n")
	cat("---------------\n")
	cat("call to SAVE:\n")
	print(object@call)
	cat("---------------\n")
	cat("Maximum likelihood estimates:\n")
	print(object@mle)
	if(!object@stage2){cat("---------------\n Message:bayesfit has not been run.\n")}
	if(object@stage2){
		cat("---------------\n")
		cat("call to bayesfit:\n")		
		print(object@bayesfitcall)
		cat("---------------\n")
		cat("Summaries of the MCMC:\n")
		print(object@mcmcsum)
	}
	cat("\n\n\n")
}

setMethod("show","summary.SAVE",
	function(object){show.summary.SAVE (object) }
)