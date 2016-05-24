###################################################################################
#' Meta Model Interface: Multi-layer Perceptron
#' 
#' Meta model based on monmlp package for multi layer perceptrons
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
#' 	spotConfig$seq.modelFit fit of the earth model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
###################################################################################
spotPredictMLP<- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"as.matrix",spotConfig$alg.roi,"monmlp","spotPredictMLP",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		xNames <- row.names(spotConfig$alg.roi)	
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]		
		y <- rawB[yNames]
		if(is.null(spotConfig$seq.mlp.n.ensemble)){spotConfig$seq.mlp.n.ensemble=20}	
		if(is.null(spotConfig$seq.mlp.hidden1)){spotConfig$seq.mlp.hidden1=2}
		if(is.null(spotConfig$seq.mlp.iter.max)){spotConfig$seq.mlp.iter.max=5000}
		if(is.null(spotConfig$seq.mlp.bag)){spotConfig$seq.mlp.bag=TRUE}
		if(is.null(spotConfig$seq.mlp.n.trials)){spotConfig$seq.mlp.n.trials=3}		
		if(is.null(spotConfig$seq.mlp.hidden2)){spotConfig$seq.mlp.hidden2=1}

		if(length(yNames)==1){
			y <- rawB[[yNames]]
			fit <- monmlp::monmlp.fit(x=as.matrix(x),y=as.matrix(y),
					hidden1=spotConfig$seq.mlp.hidden1, 
					hidden2=spotConfig$seq.mlp.hidden2, 
					iter.max=spotConfig$seq.mlp.iter.max,
					n.trials=spotConfig$seq.mlp.n.trials,
					n.ensemble=spotConfig$seq.mlp.n.ensemble,
					bag=spotConfig$seq.mlp.bag,
					silent=TRUE)	
		}
		else{#Distinction for multi criteria spot 
			y <- rawB[yNames]
			fit=list()
			res=list()
			for (i in 1:length(yNames)){
				fit[[i]] <- monmlp::monmlp.fit(x=as.matrix(x),y=as.matrix(y[,i]),
						hidden1=spotConfig$seq.mlp.hidden1, 
						hidden2=spotConfig$seq.mlp.hidden2, 
						iter.max=spotConfig$seq.mlp.iter.max,
						n.trials=spotConfig$seq.mlp.n.trials,
						n.ensemble=spotConfig$seq.mlp.n.ensemble,
						bag=spotConfig$seq.mlp.bag,
						#iter.stopped=spotConfig$seq.mlp.iter.stopped,
						silent=TRUE)
			}			
		}		
	}else{ # use existing model
		fit <-fit
	}	
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){ 	
		if(length(spotConfig$alg.resultColumn)>1){
			res=list()
			for (i in 1:length(fit)){
				res[[i]] <- monmlp::monmlp.predict(x=design,weights=fit[[i]])
			}	
		}else{res <- monmlp::monmlp.predict(x=design,weights=fit)}	
	}else{res <- NULL}	
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictMLP finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
