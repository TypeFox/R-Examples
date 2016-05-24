predict.LeafyDSA<-function(object, x, y, wt=rep(1, nrow(x)), ...){
        results <- object
	if(is.factor(y)){
		decision.rules<-results[[9]]
		}
	else{
		decision.rules<-results[[7]]
		}
	
	########HERE
	##Discuss with Karen 
	#x<-impute.test(x=x,y=y,x.test = x, y.test = y, missing = missing)
	
	#gets the max value per decision rule and saves it
	get.max.value.from.decision.rule<- function(current.decision.rule){
			length(current.decision.rule[[5]])
		}
	max.growth<-lapply(decision.rules,get.max.value.from.decision.rule)
	
	#does the prediction for partition size 1 to max value
	all.predicted.values<-lapply(decision.rules,predict, x)
			
	if(is.factor(y)){
  		y.original <- y
  		y<- ConvertFactorsToNumeric(y.original)
  			
  		#gets the predicted value based upon the max partition for that tree
  		predicted.values.for.max.growth<-mapply(function(x1,x2)  x1[[x2]]  ,all.predicted.values,max.growth, 										SIMPLIFY=FALSE)
  			  					
		convert.from.factor<-function(x){
								as.numeric(levels(x))[as.integer(x)]
							}
								
		#converts the predicted values into numbers so the computations below work
		predicted.values.for.max.growth<-lapply(predicted.values.for.max.growth,convert.from.factor)
	
		overall.pred.values <- apply(array(unlist(predicted.values.for.max.growth), dim = c(length							   (predicted.values.for.max.growth[[1]]),1, length		      				   (predicted.values.for.max.growth))), 1:2,							   MatrixMode)
	
			
		error<-categorical.error(overall.pred.values,y)
			
		#converts back to factor form so the confusion matrix shows the original type of y values
		overall.pred.values<- as.factor(Mapper(overall.pred.values,levels(y.original)))
		confusion.matrix<-create.confusion.matrix(y.original,overall.pred.values)
	
		categorical.results<- list(list("Prediction Error Rate", error),
							  list("Predicted Values", overall.pred.values),
							  list("Confusion Matrix", confusion.matrix))
		class(categorical.results)<-c('LeafyPredictions')
		return(categorical.results)

	}
	
	else{
		#gets the predicted value based upon the max partition for that tree
		predicted.values.for.max.growth<-mapply(function(x1,x2)x1												[,x2],all.predicted.values,max.growth,SIMPLIFY=FALSE)
		
		#puts in the correct for to send to numerical.predicted.value.and.error
		predicted.values.for.max.growth<-lapply(predicted.values.for.max.growth, as.matrix)
		  	
		
		calculated.results<- numerical.predicted.value.and.error(predicted.values.for.max.growth,y, wt)
		numerical.results<-list(list("Prediction Error Rate",calculated.results[[2]][[2]]),								list("Predicted Values", calculated.results[[1]][[2]]))
		class(numerical.results)<- c('LeafyPredictions')
		return(numerical.results)

		}
	
}
