categorical.predictions<- function(predicted.values.by.tree,predicted.test.set.values.by.tree,y,y.test,x,x.test,
								   y.original, y.test.original, predicted.values.by.tree.permuted){
	
	overall.pred.values <- apply(array(unlist(predicted.values.by.tree), dim = c(dim								   (predicted.values.by.tree[[1]]), length(predicted.values.by.tree))), 1:2,								   MatrixMode)
	training.set.error<- categorical.error(overall.pred.values,y)
	
	overall.pred.values<- as.factor(Mapper(overall.pred.values,levels(y.original)))
	confusion.matrix.training.set<-create.confusion.matrix(y.original,overall.pred.values)
			
			
	#if test set and training set are identical, return NA for the test set.
	test.set.error <- NA
	confusion.matrix.test.set<- NA
	overall.test.set.pred.values<-NA	
	
	if(!identical(x,x.test)){
		overall.test.set.pred.values<-apply(array(unlist(predicted.test.set.values.by.tree), dim = c(dim				    				  (predicted.test.set.values.by.tree[[1]]), 
									  length(predicted.test.set.values.by.tree))), 1:2,MatrixMode)
		test.set.error<-categorical.error(overall.test.set.pred.values,y.test)
		overall.test.set.pred.values<-as.factor(Mapper(overall.test.set.pred.values,levels(y.test.original)))
		confusion.matrix.test.set <- create.confusion.matrix(y.test.original, overall.test.set.pred.values)
	}
	
	categorical.results<- list(list("Traing Set Error:", training.set.error),
  						  list("Predicted Training Set Values:", overall.pred.values),
  						  list ("Predicted Test Set Values", overall.test.set.pred.values), 
  						  list("Test Set Error:", test.set.error),
  						  list("Training Set Confusion Matrix", confusion.matrix.training.set),
  						  list("Test Set Confusion Matrix", confusion.matrix.test.set),
  						  list("Overall Training Set Predicted Values", overall.pred.values),
  						  list("Overall Test Set Predicted Values", overall.test.set.pred.values))
  	return(categorical.results)
}

categorical.error<-function(overall.pred.values,y){
	temp <- (y != overall.pred.values)
	values.not.equal<- apply(temp, 2, sum, na.rm = TRUE)
	value.counts<- !is.na(overall.pred.values)
	total.values<- apply(value.counts,2, sum, na.rm = TRUE)
	error<- values.not.equal / total.values						   
	return(error)
	
	}

numerical.predictions<- function(predicted.values.by.tree,predicted.test.set.values.by.tree,y,y.test,x,x.test,wt,wt.test,predicted.values.by.tree.permuted){
	
	training.set.results<-numerical.predicted.value.and.error(predicted.values.by.tree,y,wt)
        #Breiman
	p <- ncol(predicted.values.by.tree.permuted[[1]])
        ntree <- length(predicted.values.by.tree.permuted)
        predicted.values.by.tree.permuted.i <- vector("list",ntree)
        training.error.permuted.difference <- rep(NA,p)
        for(i in 1:p){
          for(j in 1:ntree){
            predicted.values.by.tree.permuted.i[[j]] <- matrix(predicted.values.by.tree.permuted[[j]][,i])
          }
          results.i <- numerical.predicted.value.and.error(predicted.values.by.tree.permuted.i,y,wt)
          training.error.permuted.difference[i] <- results.i[[2]][[2]]-training.set.results[[2]][[2]]
        }
        breiman.importance.rank <- rank(-1*training.error.permuted.difference,ties.method="min")
        
	#if test set and training set are identical, return NA for the test set.        
	test.set.error<-NA #Default value if test and training set identical
	test.set.overall.pred.values<-NA #Default value if test and training set identical
	
	if(!identical(x,x.test)){
	     test.set.results<-numerical.predicted.value.and.error(predicted.test.set.values.by.tree,y.test, wt.test)
	     Predicted.Test.Set.Values<-test.set.results[[1]][[2]]
	     Test.Set.Error<-test.set.results[[2]][[2]]
	     }else{
		 Predicted.Test.Set.Values<-NA
		 Test.Set.Error<-NA     	
	}
	numerical.results<-list(list("Traininng Set Error:", training.set.results[[2]][[2]]),
  						  list("Predicted Training Set Values:", training.set.results[[1]][[2]]),
  						  list("Predicted Test Set Values",Predicted.Test.Set.Values),
  						  list("Test Set Error:",Test.Set.Error),
                                list("Training Error Permuted Difference",training.error.permuted.difference),
                                list("Breiman Importance Rank",breiman.importance.rank))
	return(numerical.results)
	}
	
numerical.predicted.value.and.error<-function(predicted.values.by.tree, y, wt){
	#Collapses all the trees into one giant sum, correctly handling the NA's
    SumResults<-Reduce('MatrixAdder',predicted.values.by.tree,right=FALSE)
    
    #The elements in the first matrix need to be 1 or NA for the MatrixCounter to run
    #correctly since it always adds 1 with each sucessive matrix. If the first matrix
    #has the actual results, the count does not start from one.
    predicted.values.by.tree[[1]]<-predicted.values.by.tree[[1]]/predicted.values.by.tree[[1]] 
    TotalCounts<-Reduce('MatrixCounter',predicted.values.by.tree,right = FALSE)
    overall.pred.values <- SumResults/TotalCounts
    
    TotalCounts<-TotalCounts/TotalCounts ## --Now we have 1 if that obs was used, NA otherwise.
	TotalCounts<-wt*TotalCounts
    tmp<- wt* (y - overall.pred.values)^2
	error<-apply(tmp, 2, sum,na.rm=TRUE) / apply(TotalCounts,2,sum,na.rm=TRUE)

	results<- list(list("predicted values", overall.pred.values),
				  list("error",error))
	return(results)
	
	}


Mapper<- function(x, key){
	mapped.predicted.values<- matrix(NA, nrow = length(x),ncol =1)
	for(i in 1:length(x)){
			#current.value<- as.numeric(levels(x[i]))[as.integer(x[i])]
			current.value<-x[i,]
			mapped.predicted.values[i,]<- key[current.value + 1]
	}
	mapped.predicted.values
}

ConvertFactorsToNumeric<- function (input){
	factor.as.numeric<-as.numeric(input)
	factor.as.numeric.adjusted<-factor.as.numeric - 1
	factor.adjusted<- as.factor(factor.as.numeric.adjusted)
	factor.adjusted
	}

create.confusion.matrix<- function(answers, predicted.values){
	
	factor.answers <- answers
	factor.predicted <- factor(predicted.values, levels = levels(factor.answers))
	confusion.matrix<-table(factor.answers, factor.predicted, useNA = "ifany")
	
	class.error<- diag(1- prop.table(confusion.matrix,1))
	confusion.matrix<-cbind(confusion.matrix,class.error)
	
	confusion.matrix

}


is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

MatrixMode <- function(x){
	cell.mode<-as.numeric(names(which.max(table(x))))
	if(length(cell.mode)==0){
		cell.mode <- NA
		}
	cell.mode
}	

MatrixAdder<-function(u,v){
	na.u<-is.na(u)
	na.v<-is.na(v)	
	ifelse(na.u & na.v, NA, ifelse(na.u, 0, u)+ ifelse(na.v,0,v))	
}

MatrixCounter<-function(u,v){
	na.u<-is.na(u)
	na.v<-is.na(v)
	ifelse(na.u & na.v, NA, ifelse(na.u, 0, u)+ ifelse(na.v,0,1))
	}

parameter.list<-function(leafy, percentage.of.variables, total.variables){

	#new.var.count<- round(percentage.of.variables*total.variables)
	new.var.count<-round(sqrt(total.variables))
	new.variables<-sample(1:total.variables,new.var.count,replace=FALSE)
	return(new.variables)
}
