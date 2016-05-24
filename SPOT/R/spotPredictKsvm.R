###################################################################################
#' Meta Model Interface: Support Vector Machine
#' 
#' Meta model based on ksvm function in kernlab package, which builds a
#' support vector machine for regression.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
###################################################################################
spotPredictKsvm <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"kernlab","spotPredictKsvm",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################
	if(is.null(fit)){		
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]		
		nx <- nrow(spotConfig$alg.roi)
		opts=list(fevals=100, reltol=1e-4)	 #optim algorithm options
		if(length(yNames)==1){
			y<-rawB[[yNames]]
			dat <- data.frame(x,y)		
			fitness<-function(xx){	
				err<-try(attributes(kernlab::ksvm(y~.,data=data.frame(x,y),epsilon=(10^xx[3]),C=(10^xx[2]),kpar=list(sigma=(10^xx[1])/nx),cross=10))$cross,silent=TRUE)
				if(class(err) == "try-error"){err<-10^4}
				err
			}			
			res <- spotOptimizationInterface(par=c(0,1,-1),fn=fitness,lower=c(-4,0,-4),upper = c(2,6,0), method="optim-L-BFGS-B", control = opts)
			fit<-try( fit<-kernlab::ksvm(y~.,data=dat,epsilon=(10^res$par[3]),C=(10^res$par[2]),kpar=list(sigma=(10^res$par[1])/nx)) ,silent=TRUE)		
			if(class(fit) == "try-error"){
				fitness<-function(xx){	
					attributes(kernlab::ksvm(y~.,data=data.frame(x,y),epsilon=(10^xx[3]),C=(10^xx[2]),kpar=list(sigma=(10^xx[1])/nx),kernel="anovadot",cross=10))$cross
				}
				res <- spotOptimizationInterface(par=c(0,1,-1),fn=fitness,lower=c(-4,0,-4),upper = c(2,6,0), method="optim-L-BFGS-B", control = opts)
				fit<-kernlab::ksvm(y~.,data=dat,epsilon=(10^res$par[3]),C=(10^res$par[2]),kpar=list(sigma=(10^res$par[1])/nx),kernel="anovadot")			
			}	
		}
		else{#Distinction for multi criteria spot 			
			fit=list()
			yy<- rawB[yNames]
			spotConfig$seq.modelFit.y<-mergedB[yNames]
			for (i in 1:length(yNames)){
				y<-yy[,i]
				dat <- data.frame(x,y)		
				fitness<-function(xx){	
					err<-try(attributes(kernlab::ksvm(y~.,data=data.frame(x,y),epsilon=(10^xx[3]),C=(10^xx[2]),kpar=list(sigma=(10^xx[1])/nx),cross=10))$cross,silent=TRUE)
					if(class(err) == "try-error"){err<-10^4}
					err
				}
				res <- spotOptimizationInterface(par=c(0,1,-1),fn=fitness,lower=c(-4,0,-4),upper = c(2,6,0), method="optim-L-BFGS-B", control = opts)
				fit[[i]]<-try( fit[[i]]<-kernlab::ksvm(y~.,data=dat,epsilon=(10^res$par[3]),C=(10^res$par[2]),kpar=list(sigma=(10^res$par[1])/nx)) ,silent=TRUE)		
				if(class(fit[[i]]) == "try-error"){
					fitness<-function(xx){	
						attributes(kernlab::ksvm(y~.,data=data.frame(x,y),epsilon=(10^xx[3]),C=(10^xx[2]),kpar=list(sigma=(10^xx[1])/nx),kernel="anovadot",cross=10))$cross
					}
					res <- spotOptimizationInterface(par=c(0,1,-1),fn=fitness,lower=c(-4,0,-4),upper = c(2,6,0), method="optim-L-BFGS-B", control = opts)
					fit[[i]]<-kernlab::ksvm(y~.,data=dat,epsilon=(10^res$par[3]),C=(10^res$par[2]),kpar=list(sigma=(10^res$par[1])/nx),kernel="anovadot")			
				}	
			}			
		}		
	}else{
		fit<-fit
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 		
		nmodel <- length(spotConfig$alg.resultColumn)
		if(nmodel>1){ #do multi criteria prediction
			res=matrix(0,nrow(design),nmodel)
			y=list()
			for (i in 1:length(fit)){ #predict		
				rres = predict(fit[[i]],design)
				rres[which(is.na(rres)==TRUE)] = median(spotConfig$alg.currentResult[,spotConfig$alg.resultColumn[i]],na.rm = TRUE)
				res[,i]= rres
			}
		}else{ #do single criteria prediction
			res <- predict(fit,design)	
			res[which(is.na(res)==TRUE)] = median(spotConfig$alg.currentResult[,spotConfig$alg.resultColumn],na.rm = TRUE) 			
		}	
	}else{res <- NULL}	
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictKsvm finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
