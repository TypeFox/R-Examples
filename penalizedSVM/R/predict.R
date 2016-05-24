predict.penSVM<-function(object, newdata, newdata.labels=NULL,labels.universe=NULL, ...){
	# calculate test error, confusion table, sensitivity and specificity for test data
	# input:
	# newdata: matrx of ratios, rows - samples, columns - clones
	# newdata.labels - vector of classes for samples
	# object - fit of trained data (producedby svm.fs) 
	# labels.universe: important for models produced by loocv: all possible labels in the particular data set
	 	
	pred.class<- NULL
	tab.classes<- NULL
	error.classes<- NA
	sensitivity.classes<- NA
	specificity.classes<- NA
	
	# fit model 	
	f<-object$model	
	# if we have a model, calulate the test error...
	if (length(f)<=1 )  {
		stop("Model is empty!")
	}else {
		# 2 different data possible: for GENES and for CLASSES (here only for classes)
		
		# separating line: x*w+b=f(x)
		sep = as.matrix(newdata)[,f$xind, drop=FALSE] %*% f$w + f$b
		
		# for classes -1 and 1 :  class = 2*(sep > 0) - 1
		pred.class = factor(2*as.numeric (sep > 0) -1)
				
		if (!is.null(newdata.labels)){	
		
			# missclassification table for CLASSES		
			tab.classes<-table(pred.class, newdata.labels)
			
					
		
		# 3. sensitivity and specificity  for CLASSES
		if (!is.null(tab.classes)){
			# if in test only one sample --> extend tab.classes to complete labels to  labels.universe	
			# example : 
			#tab.classes
      #    newdata.labels
			#pred.class 1
      #   1 1
			# labels.universe= c("-1","1")
		
			if (nrow(tab.classes)!= ncol(tab.classes) | nrow(tab.classes)== 1 ) 	  tab.classes<-.extend.to.quad.matrix (tab.classes, labels.universe=labels.universe)
			# sensitivity = TP/ all P = TP /(TP + FN)
			sensitivity.classes<- tab.classes[2,2]/sum(tab.classes[,2])
			# specificity = TN/ all N = TN /(TN + FP)
			specificity.classes <- tab.classes[1,1]/sum(tab.classes[,1])
			# secondary diagonal
			sec.diag<-c(); 	for (j in 1:ncol( tab.classes)) sec.diag<-  c(sec.diag,  tab.classes[ncol( tab.classes)-j+1,j]  )
			error.classes<- ( sum(sec.diag) ) / sum( tab.classes)
		}
	}# end of  if (!is.null(newdata.labels))
	
	return(list(pred.class = pred.class,
							fitted=sep,
							tab=tab.classes,
							error=error.classes,
							sensitivity=sensitivity.classes,
							specificity=specificity.classes ))
	} # end of ifelse (is.null(f))
}
