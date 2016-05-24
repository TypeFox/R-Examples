AbsoluteQuantification <- function(data, ...) UseMethod("AbsoluteQuantification")

AbsoluteQuantification.default <- function(data, total_protein_concentration=1, ...) {
    if (length(setdiff(c("run_id","protein_id","concentration","response"),names(data))) > 0) stop("No suitable input for AbsoluteQuantification. Did you use ProteinInference on your data.frame?")	

	object = list()

	data.selection<-data
	
	data.selection$normalized_response <- data.selection$response / sum(data.selection$response)
	data.selection$response <- log(data.selection$response)
	data.selection$normalized_response <- log(data.selection$normalized_response)
	
	object$is_calibrated <- FALSE
	
	if (dim(subset(data.selection,data.selection$concentration != "?"))[1] > 2) {	
		object$is_calibrated <- TRUE
		object$calibration <- subset(data.selection,data.selection$concentration != "?")
		object$calibration$concentration <- log(as.numeric(object$calibration$concentration))
		object$prediction <- data.selection
		object$prediction$concentration<-"?"
				
		object$model <- lm(concentration ~ response,object$calibration)
		object$mfe <- mean(folderror.AbsoluteQuantification(exp(predict(object$model,object$calibration)),exp(object$calibration$concentration)))
		object$r.squared <- summary(object$model)$r.squared
		object$calibration_covar <- (100*sd(object$calibration$response)) / mean(object$calibration$response)
	}
	
	data.selection$normalized_concentration <- data.selection$normalized_response + log(total_protein_concentration)
	object$estimation <- data.selection[,c("run_id","protein_id","normalized_response","normalized_concentration")]
	object$normalized_concentration_covar <- (100*sd(data.selection$normalized_response)) / mean(data.selection$normalized_response)

	class(object) <- "AbsoluteQuantification"
	return(object)
}

predict.AbsoluteQuantification <- function(object, ...) {
    if (!inherits(object, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if (!object$is_calibrated) stop("Method not supported for AbsoluteQuantification objects without calibration proteins.")
	object$prediction$concentration <- predict(object$model,newdata=object$prediction)
	return(object)
}

folderror.AbsoluteQuantification <- function(predicted,true, ...) {
    fea = rep(0,0)
    for(i in 1:length(predicted)) {
    	fe <- predicted[i] / true[i]
		if (fe < 1) fe = 1 / fe
		fea = c(fea, fe) 
	}
	return(fea)
}

cval <- function(object, ...)  UseMethod("cval")

cval.default <- function(object, ...)
    stop("No method implemented for this class of object")

cval.AbsoluteQuantification <- function(object,cval_method="mc",mcx=1000, ...) {
    if (!inherits(object, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if (!object$is_calibrated) stop("Method not supported for AbsoluteQuantification objects without calibration proteins.")
	rs_sample <- 0
	mfe_sample <- 0
	object$cval$rs_array <- c()
	object$cval$mfe_array <- c()
	j <- 1
	
	if (cval_method=="mc") {
		while (j <= mcx) {
			partition1 <- sample(c(1:dim(object$calibration)[1]),floor((1/3)*dim(object$calibration)[1]))
			partition2 <- sample(c(1:dim(object$calibration)[1])[-partition1],floor((1/3)*dim(object$calibration)[1]))
			partition3 <- c(1:dim(object$calibration)[1])[-c(partition1,partition2)]
			
			regression1 <- lm(formula(concentration ~ response),object$calibration[c(partition1,partition2),])
			rs1 <- 1 - sum((abs(object$calibration[partition3,]$concentration - predict(regression1,object$calibration[partition3,])))^2) / sum((abs(object$calibration[partition3,]$concentration - mean(object$calibration[partition3,]$concentration)))^2)
			mfe1 <- mean(folderror.AbsoluteQuantification(as.vector(exp(predict(regression1,object$calibration[partition3,]))),exp(object$calibration[partition3,]$concentration)))
			regression2 <- lm(formula(concentration ~ response),object$calibration[c(partition1,partition3),])
			rs2 <- 1 - sum((abs(object$calibration[partition2,]$concentration - predict(regression2,object$calibration[partition2,])))^2) / sum((abs(object$calibration[partition2,]$concentration - mean(object$calibration[partition2,]$concentration)))^2)
			mfe2 <- mean(folderror.AbsoluteQuantification(as.vector(exp(predict(regression2,object$calibration[partition2,]))),exp(object$calibration[partition2,]$concentration)))
			regression3 <- lm(formula(concentration ~ response),object$calibration[c(partition2,partition3),])
			rs3 <- 1 - sum((abs(object$calibration[partition1,]$concentration - predict(regression3,object$calibration[partition1,])))^2) / sum((abs(object$calibration[partition1,]$concentration - mean(object$calibration[partition1,]$concentration)))^2)
			mfe3 <- mean(folderror.AbsoluteQuantification(as.vector(exp(predict(regression3,object$calibration[partition1,]))),exp(object$calibration[partition1,]$concentration)))
			
			object$cval$mfe_array <- append(object$cval$mfe_array,c(mfe1,mfe2,mfe3))
			object$cval$rs_array <- append(object$cval$rs_array,c(rs1,rs2,rs3))
			
			j <- j + 1
		}
		object$cval$r.squared <- mean(object$cval$rs_array[which(is.finite(object$cval$rs_array))])
		object$cval$mfe <- mean(object$cval$mfe_array[which(is.finite(object$cval$mfe_array))])
	}
	else if (cval_method=="boot") {
		while (j <= mcx) {
			training <- sample(dim(object$calibration)[1],replace=TRUE)
			
			regression <- lm(formula(concentration ~ response),object$calibration[training,])
			rs <- 1 - sum((abs(object$calibration[-unique(training),]$concentration - predict(regression,object$calibration[-unique(training),])))^2) / sum((abs(object$calibration[-unique(training),]$concentration - mean(object$calibration[-unique(training),]$concentration)))^2)
			fe <- folderror.AbsoluteQuantification(as.vector(exp(predict(regression,object$calibration[-training,]))),exp(object$calibration[-training,]$concentration))

			object$cval$mfe_array <- append(object$cval$mfe_array,mean(fe))
			object$cval$rs_array <- append(object$cval$rs_array,rs)
			
			j <- j + 1
		}
		object$cval$r.squared <- mean(object$cval$rs_array[which(is.finite(object$cval$rs_array))])
		object$cval$mfe <- mean(object$cval$mfe_array[which(is.finite(object$cval$mfe_array))])
	}
	else if (cval_method=="loo") {
		j <- 1
		while (j <= dim(object$calibration)[1]) {
			training <- c(1:dim(object$calibration)[1])[-j]
			test <- j

			regression1 <- lm(formula(concentration ~ response),object$calibration[training,])
			mfe1 <- 1 + abs(1 - abs(exp(predict(regression1,object$calibration[test,])) - exp(object$calibration[test,]$concentration)) / exp(object$calibration[test,]$concentration))
			
			mfe_sample <- mfe_sample + mfe1

			j <- j + 1
		}
		object$cval$r.squared <- summary(object$model)$r.squared
		object$cval$mfe <- mfe_sample / dim(object$calibration)[1]
	}
	return(object)
}

print.AbsoluteQuantification <- function(x, ...) {
    if (!inherits(x, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
	cat("AbsoluteQuantification\n")
	cat("\n")
	cat("Number of proteins: ")
	cat(length(unique(x$estimation$protein_id)))
	cat("\n")
	cat("Calibration Trainingset size: ")
	cat(dim(x$calibration)[[1]])
	cat("\n")
	cat("Calibration Testset size: ")
	cat(dim(x$prediction)[[1]])
	cat("\n")
	cat("Calibration Regression mean-fold error: ")
	cat(x$mfe)
	cat("\n")
	cat("Calibration Regression R-squared: ")
	cat(x$r.squared)
	cat("\n")
	cat("Calibration cross-validation mean-fold error: ")
	cat(x$cval$mfe)
	cat("\n")
	cat("Calibration cross-validation R-squared: ")
	cat(x$cval$r.squared)
	cat("\n")
	cat("Calibration CV: ")
	cat(x$calibration_covar)
	cat("\n")
	cat("Normalized Concentration CV: ")
	cat(x$normalized_concentration_covar)
	cat("\n")
}

plot.AbsoluteQuantification <- function(x, ...) {
	concentration <- response <- NULL

    if (!inherits(x, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if (!x$is_calibrated) stop("Method not supported for AbsoluteQuantification objects without calibration proteins.")
	if (dim(subset(x$prediction,concentration == "?"))[1] > 0) {
		plot(data.frame(log10(exp(x$calibration$response)),log10(exp(as.numeric(x$calibration$concentration)))),col="red",xlab='log10(intensity)',ylab='log10(concentration)',xlim=c(min(log10(exp(as.numeric(x$calibration$response)))),max(log10(exp(as.numeric(x$calibration$response))))),ylim=c(min(log10(exp(as.numeric(x$calibration$concentration)))),max(log10(exp(as.numeric(x$calibration$concentration))))))
		log10_model <-x$model
		log10_model$coefficients[1] <- log10(exp(log10_model$coefficients[1]))
		abline(log10_model, col="red")
		if (is.null(x$cval)) {
			title(paste('MFE:',format(x$mfe, digits=4),'RSQ:',format(x$r.squared, digits=4)))
		}
		else {
			title(paste('CV-MFE:',format(x$cv$mfe, digits=4),'CV-RSQ:',format(x$cv$r.squared, digits=4)))
		}
	}
	else if (dim(subset(x$prediction,concentration == "?"))[1] == 0) {
		calibration.df <- data.frame(x$calibration$response,as.numeric(x$calibration$concentration))
		names(calibration.df) <- c("response","concentration")
		prediction.df <- data.frame(x$prediction$response,as.numeric(x$prediction$concentration))
		names(prediction.df) <- c("response","concentration")
		merged <- rbind(calibration.df,prediction.df)
		merged <- subset(merged,is.finite(concentration) & is.finite(response))
		plot(log10(exp(prediction.df)),xlab='log10(intensity)',ylab='log10(concentration)',xlim=c(min(log10(exp(as.numeric(merged$response)))),max(log10(exp(as.numeric(merged$response))))),ylim=c(min(log10(exp(as.numeric(merged$concentration)))),max(log10(exp(as.numeric(merged$concentration))))))
		points(log10(exp(calibration.df)), col ="red")
		log10_model <-x$model
		log10_model$coefficients[1] <- log10(exp(log10_model$coefficients[1]))
		abline(log10_model, col="red")
		if (is.null(x$cval)) {
			title(paste('MFE:',format(x$mfe, digits=4),'RSQ:',format(x$r.squared, digits=4)))
		}
		else {
			title(paste('CV-MFE:',format(x$cv$mfe, digits=4),'CV-RSQ:',format(x$cv$r.squared, digits=4)))
		}
	}
}

hist.AbsoluteQuantification <- function(x, ...) {
    if (!inherits(x, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if (!x$is_calibrated) stop("Method not supported for AbsoluteQuantification objects without calibration proteins.")
    if (is.null(x$cval)) stop("Apply cval method to object before plotting the histogram.")

	hist(x$cval$mfe_array, freq=FALSE,main=paste('Histogram of MFE\nMean = ',format(mean(x$cval$mfe_array), digits=4), ", 95% CI =", format(1.96*sd(x$cval$mfe_array), digits=4) ), breaks = 40, xlab="mean fold error")
	xfit = seq(min(x$cval$mfe_array), max(x$cval$mfe_array),length=40)
	yfit = dnorm(xfit, mean=mean(x$cval$mfe_array), sd=sd(x$cval$mfe_array))
	lines(xfit, yfit, col="red")
}

export <- function(x, file, ...)  UseMethod("export")

export.default <- function(x, file, ...)
    stop("No method implemented for this class of object")

export.AbsoluteQuantification <- function(x, file, ...) {
    if (!inherits(x, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if (x$is_calibrated) {
    	data <- merge(x$prediction[,c("protein_id","response","concentration")],x$estimation, all=TRUE)
   	    data$response <- exp(data$response)
  		data$concentration <- exp(as.numeric(data$concentration))
    }
    else {
    	data <- x$estimation
    }
            
    data$normalized_response <- exp(data$normalized_response)
    data$normalized_concentration <- exp(data$normalized_concentration)
    
	write.csv(data, file = file, ...)
}

pivot <- function(x, ...)  UseMethod("pivot")

pivot.default <- function(x, ...) {
    if (!inherits(x, "data.frame")) stop("Is not a data.frame")

	if ("sec_id" %in% names(x)) {
		data<-acast(x, protein_id ~ sec_id, value.var="response", fill=0)
	}
	else {
		data<-acast(x, protein_id ~ run_id, value.var="response", fill=0)
	}	

	return(data)
}

pivot.AbsoluteQuantification <- function(x, ...) {
    if (!inherits(x, "AbsoluteQuantification")) stop("Is not a AbsoluteQuantification object")
    if ("?" %in% x$prediction$concentration) stop("Apply predict before pivot to AbsoluteQuantification object")

    if (x$is_calibrated) {
    	x$prediction$concentration<-exp(as.numeric(x$prediction$concentration))
		tdata<-acast(x$prediction, protein_id ~ run_id, value.var="concentration", fill=0)
    }
    else {
     	x$prediction$normalized_concentration<-exp(x$prediction$normalized_concentration)
		tdata<-acast(x$estimation, protein_id ~ run_id, value.var="normalized_concentration", fill=0)
    }

	return(tdata)
}
