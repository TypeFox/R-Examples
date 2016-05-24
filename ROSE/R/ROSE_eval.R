#Last modified on 01/30/2014 

ROSE.eval <- function(formula, data, learner, acc.measure="auc", extr.pred=NULL, method.assess="holdout", K=1, B=100, control.rose=list(), control.learner=list(), control.predict=list(), control.accuracy=list(), trace=FALSE, subset=options("subset")$subset, na.action=options("na.action")$na.action, seed) 
{
		#check arguments: formula and learner are mandatory 
		if( missing(formula) ) 
			stop("A formula is required.\n")
		if( missing( learner ) ) 
			stop("Argument 'learner' is missing, with no default. \n")

		#check if provided learner is "standard" in the sense that it has an associated predict method with arguments "object" and "newdata"
		func.name <- as.character(substitute(learner))
		if( any( methods(class=func.name)==paste("predict.",func.name,sep="") ) )
			flg.learner <- 1
		else
			flg.learner <- 0

	mc <- match.call()
	formula.env <- attr(formula,".Environment")

	varnames.func <- all.vars(formula, functions=TRUE)
	varnames <- all.vars(formula, functions=FALSE)

###catch the original data.frame/variables in formula.env
	name.data <- NULL

		if(missing(data))
		{
				if( any( varnames.func%in%c("$","[","]")) )
					name.data <- varnames[1]
		}
		else		
			name.data <- as.character(mc$data)

	
		##this is the case for variables not contained in a data frame in formula.env
		if(is.null(name.data))
			name.data <- varnames

	###store original data.frame/variables in formula.env
	data.orig <- sapply(name.data, function(x) get(x, envir=formula.env) )
	cn.order.orig <- attributes(data.orig)$dimnames[[1]]
	###end

	###keep formula unchanged for the learner
	formula.learn <- formula
	###drop trasformations etc from formula to provide a nice formula to ROSE
	formula.rose <- adj.formula(formula, data)

		if(missing(data))
			lst.model.frame <- list(formula=formula.rose, data=NULL, subset=subset, na.action=na.action)
		else
			lst.model.frame <- list(formula=formula.rose, data=data, subset=subset, na.action=na.action)

	###create data set for ROSE and prediction
#	mc$formula <- formula.rose
#	m <- match(c("formula", "data", "na.action", "subset"), names(mc), 0L)
#	mf <- mc[c(1L, m)]
#	mf[[1L]] <- as.name("model.frame")
#	mf <- eval(mf, parent.frame())
	mf <- do.call(model.frame,lst.model.frame)
	cn <- rownames( attributes( attributes(mf)$terms )$factors )
	data.st <- data.frame(mf)
	y <- data.st[,1]

		if( any( varnames.func%in%c("$")) ) 
			colnames(data.st) <- gsub(paste(name.data, ".", sep=""), "", colnames(data.st))

	##create new formula for ROSE with the right environment
	formula.rose <- formula(data.st)
	###end 

		#right order of columns as specified in data
	d <- NCOL(data.st)-1
		if( !missing(data) )
				if(d!=1)
					data.st <- data.st[cn.order.orig]

	#check accuracy estimator
	method.assess <- match.arg(method.assess, choices=c("holdout", "LKOCV", "BOOT"))
		if(!method.assess %in% c("holdout", "LKOCV", "BOOT") ) stop("Method for model assessment must be one among 'holdout', 'BOOT' or 'LKOCV'.\n") 

	#check accuracy measure 
	acc.measure <- match.arg(acc.measure, choices=c("auc", "precision", "recall", "F"))
		if(!acc.measure %in% c("auc", "precision", "recall", "F") ) stop("Accuracy measure must be one among 'auc' 'precision', 'recall' or 'F'.\n") 
		if(acc.measure=="auc")
		{
			fun.accuracy <- roc.curve 
			##as default, do not plot the roc curve in ROSE.eval when the user do not want it
			if( is.null(names(control.accuracy)) ) control.accuracy <- c(control.accuracy, list("plotit"=FALSE))
			if( !names(control.accuracy)=="plotit" ) control.accuracy <- c(control.accuracy, list("plotit"=FALSE))
		}
		else 
		{
			fun.accuracy <- accuracy.meas
		}

	pos.accuracy <- match(acc.measure,c("auc","precision", "recall", "F")) + 1
	method.assess.inn <- method.assess

		if(!missing(seed)) set.seed(seed)

		if(trace)
		{
			ind <- ifelse( B<50, 1, ifelse( B<500, 10, 100 ) )
			cat("Iteration:", "\n")
		}

		if( method.assess.inn =="holdout" )
		{
			method.assess.inn <- "BOOT"
			B <- 1
		}

		if( method.assess.inn=="BOOT" )
		{
			if(trace) max.ind <- floor(B/ind)*ind
			acc.vec <- vector(mode="numeric", length=B) 

			if( flg.learner )
			{
				#functions with "standard" behaviour
				for(i in 1:B)
				{
					data.rose <- do.call(ROSE, c(list(formula=formula.rose, data=data.st), control.rose))$data
					fit <- do.call(learner, c(list(formula=formula.learn, data=data.rose), control.learner))
					pred <- do.call(predict, c(list(object=fit, newdata=data.st), control.predict))
						if(!is.null(extr.pred)) pred <- extr.pred(pred)
					acc.vec[i] <- do.call(fun.accuracy, c(list(response=y, predicted=pred), control.accuracy))[[pos.accuracy]]
						if(trace) if(i %% ind == 0) {if( i!=max.ind ) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
				}
			}
			else
			{
				#user defined functions with "non-standard" behaviour
				for(i in 1:B)
				{
					data.rose <- do.call(ROSE, c(list(formula=formula.rose, data=data.st), control.rose))$data
					pred <- do.call(learner, c(list(data=data.rose, newdata=data.st[,cn[-1]]), control.learner))
					acc.vec[i] <- do.call(fun.accuracy, c(list(response=y, predicted=pred), control.accuracy))[[pos.accuracy]]
						if(trace) if(i %% ind == 0) {if( i!=max.ind ) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
				}
			}
		}
		else
		{
			pred <- y.cp <- numeric(0)

				if(trace)	max.ind <- floor(B/ind)*ind
			#n.obs to leave out
			if(K%%1!=0) stop("Leave K out CV: K must be an integer\n")
			n.g <- K

			#number of subsets
			if(length(data.st[,1])%%n.g==0)
			{
				K <- length(data.st[,1])/n.g
				ind.g <- sample( rep(1:K, n.g) )
			}
			else
			{
				K <- floor(length(data.st[,1])/n.g) + 1
				n.g.remain <- length(data.st[,1])-floor(length(data.st[,1])/n.g)*n.g
					message(paste("\nLeave K out CV: the sample size is not a multiple of K. \nThe routine has automatically created", K-1, "subsets of size", n.g, "and one subset of size", n.g.remain,"."))
				ind.g <- sample( c(rep(1:(K-1), n.g), rep(K,n.g.remain) ) )
			}

			B <- K

			if( flg.learner )
			{
				#functions with "standard" behaviour
				for(i in 1:B)
				{
					data.rose <- do.call(ROSE, c(list(formula=formula.rose, data=data.st[ -which(ind.g==i) ,]), control.rose))$data
					fit <- do.call(learner, c(list(formula=formula.learn, data=data.rose), control.learner))
					predi <- do.call(predict, c(list(object=fit, newdata=data.st[ which(ind.g==i) ,]), control.predict))
						if(!is.null(extr.pred)) predi <- extr.pred(predi)
					pred <- c(pred,predi)
						if(trace) if(i %% ind == 0) {if( i!=max.ind ) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
					y.cp <- c(y.cp,y[which(ind.g==i)])
				}
				acc.vec <- do.call(fun.accuracy, c(list(response=y.cp, predicted=pred), control.accuracy))[[pos.accuracy]]
			}
			else
			{
				#user defined functions with "non-standard" behaviour
				for(i in 1:B)
				{
					data.rose <- do.call(ROSE, c(list(formula=formula.rose, data=data.st[ -which(ind.g==i) ,]), control.rose))$data
					predi <- do.call(learner, c(list(data=data.rose, newdata=data.st[which(ind.g==i),cn[-1]]), control.learner))
					pred <- c(pred,predi)
						if(trace) if(i %% ind == 0) {if( i!=max.ind ) cat(i, ", ", sep="") else cat(i, "\n", sep="")}  
					y.cp <- c(y.cp,y[which(ind.g==i)])
				}
				acc.vec <- do.call(fun.accuracy, c(list(response=y.cp, predicted=pred), control.accuracy))[[pos.accuracy]]
			}
		}

#		out <- list(Call = match.call(), method=method.assess, measure = acc.measure, acc = acc.vec)
		out <- list(Call = mc, method=method.assess, measure = acc.measure, acc = acc.vec)
		class(out) <- "ROSE.eval"
		out
}

##print method for ROSE.eval
print.ROSE.eval <- function(x, ...) 
{
		if (x$method =="BOOT")		method <- "Bootstrap"
		if (x$method =="LKOCV")		method <- "Leave K out cross-validation"
		if (x$method =="holdout")	method <- "Holdout"

	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")

		if (method == "Bootstrap")
			cat( paste(method, " estimate of ", x$measure, " on ", length(x$acc), " samples: \n ", sep="") )
		else
			cat( paste(method, " estimate of ", x$measure, ": ", sep="") )

	cat(sprintf("%.3f",x$acc),"\n")
}

###summary method for ROSE.eval
summary.ROSE.eval <- function(object, ...) 
{
	acc<-object$acc
		if (length(acc) > 1) acc <- summary(acc) 
	LST <- list( call=object$Call, method=object$method, measure=object$measure, acc=acc ) 
	class(LST) <- "summary.ROSE.eval"
	LST
}

###print method for summary
print.summary.ROSE.eval <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$call)
	cat("\n")

		if (x$method =="BOOT") method <- "Bootstrap"
		if (x$method =="LKOCV") method <- "Leave K out cross-validation"
		if (x$method =="holdout") method <- "Holdout"

		if(x$method !="BOOT")
			cat( paste(method, " estimate of ", x$measure, ": ", sprintf("%.3f",x$acc),"\n", sep="") )
		else
		{
			cat( "Summary of bootstrap distribution of auc: \n" )
			print(x$acc)
			cat("\n")
		}
}
