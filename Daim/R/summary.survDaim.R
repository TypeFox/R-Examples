

DaimList <- function(x,times){
	x <- data.frame(x)
	Names <- row.names(x)
	Names[1] <- "mean"
	ans <- data.frame(times,t(x),row.names=NULL)
	names(ans) <- c("time",Names)
	ans
}




summary.survDaim <- function(object, sd=FALSE, CI.prob = c(0.025,0.975), ...)
{
	meth <- unlist(object$method)
	meth.names <- names(meth)
	meth <- as.vector(meth)
	times <- object$times
	if(class(object)[2] != "cv"){
		if(!sd){
			erg.iauc <- lapply(object[4:7], function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.632p <- lapply(object$auc.632p, function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.632 <- lapply(object$auc.632, function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.oob <- lapply(object$auc.oob, function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.app <- lapply(object$auc.app, function(x) c(mean(x),quantile(x,prob=CI.prob)))
		}else{
			erg.iauc <- lapply(object[4:7], function(x) c(mean(x),sd(x)))
			erg.auc.times.632p <- lapply(object$auc.632p, function(x) c(mean(x),sd(x)))
			erg.auc.times.632 <- lapply(object$auc.632, function(x) c(mean(x),sd(x)))
			erg.auc.times.oob <- lapply(object$auc.oob, function(x) c(mean(x),sd(x)))
			erg.auc.times.app <- lapply(object$auc.app, function(x) c(mean(x),sd(x)))
		}
		erg.iauc <- data.frame(erg.iauc)
		erg.auc.times <- list(auc.times.632p=DaimList(erg.auc.times.632p,times), auc.times.632=DaimList(erg.auc.times.632,times), 
							  auc.times.oob=DaimList(erg.auc.times.oob,times), auc.times.app=DaimList(erg.auc.times.app,times))
		row.names(erg.iauc)[1] <- "mean"
		
		ans <- list(call=object$formula,
					method=list(meth.names=meth.names,meth=meth),
					IAUC=erg.iauc,
					AUC.times = erg.auc.times,
					time=times)
	}
	else{
		if(!sd){
			erg.iauc <- lapply(object[6:7], function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.cv <- lapply(object$auc.oob, function(x) c(mean(x),quantile(x,prob=CI.prob)))
			erg.auc.times.app <- lapply(object$auc.app, function(x) c(mean(x),quantile(x,prob=CI.prob)))
		}else{
			erg.iauc <- lapply(object[6:7], function(x) c(mean(x),sd(x)))
			erg.auc.times.cv <- lapply(object$auc.oob, function(x) c(mean(x),sd(x)))
			erg.auc.times.app <- lapply(object$auc.app, function(x) c(mean(x),sd(x)))
		}
		erg.iauc <- data.frame(erg.iauc)
		erg.auc.times <- list(auc.times.cv=DaimList(erg.auc.times.cv,times), auc.times.app=DaimList(erg.auc.times.app,times))
		row.names(erg.iauc)[1] <- "mean"
		
		ans <- list(call=object$formula,
					method=list(meth.names=meth.names,meth=meth),
					IAUC=erg.iauc,
					AUC.times = erg.auc.times,
					time=times)
	}
	class(ans) <- c("summary.survDaim",class(object)[2])
	ans
}


print.summary.survDaim <- function(x, digits=max(3, getOption("digits") - 3), ...)
{
	method <- paste(paste(x$method$meth.names,x$method$meth,sep=" = "),"",sep=",")
	method[length(method)] <- gsub(",",".",method[length(method)])
	cat("\nPerformance obtained by:\n")
	cat("\nCall:\n ")
	print(formula(x$call))
	cat("\nsurvDaim parameters: \n")
	cat(method,labels=" ",fill=TRUE)
	cat("\nResult: \n")
	IAUC <- format(round(x$IAUC,digits=digits),digits=digits)	
	erg.iauc.mean <- format(x$IAUC[1,],digits=digits)
	erg.iauc.CI <- x$IAUC[-1,]
	erg.iauc.CI.names <- row.names(erg.iauc.CI)
	
	erg.iauc.CI <- paste("(",apply(erg.iauc.CI,2,
			function(x) paste(format(round(as.numeric(x),digits),digits=digits),sep="",collapse=", ")),")",sep="")
	erg <- data.frame(mean=as.numeric(erg.iauc.mean), CI=erg.iauc.CI)
	names(erg) <- c("mean", paste("(",paste(erg.iauc.CI.names,collapse=", "),")   ",sep=""))
	if(class(x)[2] != "cv"){
		row.names(erg) <- c(" 0.632+    "," 0.632"," oob"," app")

		best.AUC.times.632p <- which.max(x$AUC.times$auc.times.632p[,2])
		best.AUC.times.632 <- which.max(x$AUC.times$auc.times.632[,2])
		best.AUC.times.oob <- which.max(x$AUC.times$auc.times.oob[,2])
		best.AUC.times.app <- which.max(x$AUC.times$auc.times.app[,2])
		cat("------------------------------------\n")
		cat(" Int. AUC                           \n")
		cat("------------------------------------\n")
		print(format(erg,digits=digits))
		cat("====================================\n")
		cat("\n\nBest time-dependent AUC:\n",sep="")
		cat("------------------------------------------------\n")
		cat("Est.method |", paste(names(x$AUC.times$auc.times.oob)," ",sep=" | "), "\n")
		cat("------------------------------------------------\n")
		
		tmp <- as.numeric(x$AUC.times$auc.times.632p[best.AUC.times.632p,])
		cat(" 0.632+    |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		tmp <- as.numeric(x$AUC.times$auc.times.632[best.AUC.times.632,])
		cat(" 0.632     |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		tmp <- as.numeric(x$AUC.times$auc.times.oob[best.AUC.times.oob,])
		cat(" oob       |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		tmp <- as.numeric(x$AUC.times$auc.times.app[best.AUC.times.app,])
		cat(" app       |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		cat("------------------------------------------------\n")
		cat("\n")
		invisible(x)
	}
	else{
		row.names(erg) <- c(" oob"," app")

		best.AUC.times.cv <- which.max(x$AUC.times$auc.times.cv[,2])
		best.AUC.times.app <- which.max(x$AUC.times$auc.times.app[,2])
		
		cat("------------------------------------\n")
		cat(" Int. AUC                           \n")
		cat("------------------------------------\n")
		print(format(erg,digits=digits))
		cat("====================================\n")
		cat("\n\nBest time-dependent AUC:\n",sep="")
		cat("------------------------------------------------\n")
		cat("Est.method |", paste(names(x$AUC.times$auc.times.cv)," ",sep=" | "), "\n")
		cat("------------------------------------------------\n")
		
		tmp <- as.numeric(x$AUC.times$auc.times.cv[best.AUC.times.cv,])
		cat(" oob       |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		tmp <- as.numeric(x$AUC.times$auc.times.app[best.AUC.times.app,])
		cat(" app       |",paste(c(tmp[1],format(tmp[-1],digits=digits,nsmall=digits))," ",sep=" |"), "\n")
		cat("------------------------------------------------\n")
		cat("\n")
		invisible(x)
	}
}

