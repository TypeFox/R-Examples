APEX <- function(data, ...) UseMethod("APEX")

APEX.default <- function(data, ...) {
	class(data) <- "data.frame"

	object = list()
	object$classes <- all.vars(formula)
	object$training <- data
				
	object$model <- randomForest(apex ~ ., data=object$training[,-which(colnames(object$training) == "peptide_sequence")], classwt=table(object$training$apex)/length(object$training$apex))
	
	class(object) <- "APEX"
	
	return(object)
}

predict.APEX <- function(object, newdata=NULL, ...) {
    if (!inherits(object, "APEX")) stop("Is not a APEX object")
	class(newdata) <- "data.frame"
	object$prediction <- newdata
	object$prediction$apex <- predict(object$model,object$prediction,type="prob")[,2]
					
	return(object)
}

cval.APEX <- function(object, folds=10, ...) {	
    if (!inherits(object, "APEX")) stop("Is not a APEX object")
	object$cval <- list()
	
	folds.list <- createFolds(array(1:dim(object$training)[1]), k = folds, list = TRUE)
	
	i <- 1
	while (i <= folds) {
		model <- randomForest(apex ~ ., data=object$training[-folds.list[[i]],-which(colnames(object$training) == "peptide_sequence")], classwt=table(object$training$apex)/length(object$training$apex))
		object$cval$predicted[folds.list[[i]]] <- predict(model,object$training[folds.list[[i]],],type="prob")[,2]

		i <- i+1
	}
	object$cval$apex <- object$training$apex
	object$cval$roc.prediction <- prediction(object$cval$predicted,object$cval$apex)
	object$cval$roc.performance <- performance(object$cval$roc.prediction,"tpr","fpr")
	object$cval$roc.auc <- performance(object$cval$roc.prediction,"auc")
	
	return(object)
}

print.APEX <- function(x, ...) {
    if (!inherits(x, "APEX")) stop("Is not a APEX object")
	cat("APEX\n")
	cat("Features: ")
	cat(dim(x$training)[[2]])
	cat("\n")
	cat("Trainingset size: ")
	cat(dim(x$training)[[1]])
	cat("\n")
	cat("Testset size: ")
	cat(dim(x$prediction)[[1]])
	cat("\n")
	if (!is.null(x$cval)){
		cat("AUC: ")
		cat(x$cval$roc.auc@y.values[[1]])
		cat("\n")
	}
}

plot.APEX <- function(x, ...) {
    if (!inherits(x, "APEX")) stop("Is not a APEX object")
	if (!is.null(x$cval)) {
		plot(x$cval$roc.performance)
	}
}
