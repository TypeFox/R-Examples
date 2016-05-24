# last modified 26 December 2010 by J. Fox

CoxZPH <- function(){
	.activeModel <- ActiveModel()
	command <- paste(".CoxZPH <- cox.zph(", .activeModel, ")", sep="")
	doItAndPrint(command)
	doItAndPrint(".CoxZPH")
	nvar <- ncol(.CoxZPH$y)
	doItAndPrint(paste(".b <- coef(", .activeModel, ")", sep=""))
	doItAndPrint(paste(".mfrow <- par(mfrow = mfrow(", nvar, "))", sep=""))
	for (i in 1:nvar){
		doItAndPrint(paste("plot(.CoxZPH[", i, "])", sep=""))
		doItAndPrint(paste("abline(h=.b[", i, "], lty=3)", sep=""))
		doItAndPrint(paste("abline(lm(.CoxZPH$y[,", i, "] ~ .CoxZPH$x), lty=4)", sep=""))
	}
	doItAndPrint("par(mfrow=.mfrow)")
	logger("remove(.CoxZPH, .mfrow, .b)")
	remove(.CoxZPH, .mfrow, .b, envir=.GlobalEnv)
}

CoxDfbetas <- function(){ # works for survreg models as well
	command <- paste(".dfbetas <- as.matrix(residuals(", ActiveModel(), ', type="dfbetas"))', sep="")
	doItAndPrint(command)
	command <- if (coxphP())
			paste("colnames(.dfbetas) <- names(coef(", ActiveModel(), "))", sep="")
		else paste("colnames(.dfbetas) <- rownames(summary(", ActiveModel(), ")$table)", sep="")
	doItAndPrint(command)
	ncol <- ncol(.dfbetas)
	doItAndPrint(paste(".mfrow <- par(mfrow = mfrow(", ncol, "))", sep=""))
	for (col in colnames(.dfbetas)){
		doItAndPrint(paste('plot(.dfbetas[,"', col, '"])', sep=""))
		doItAndPrint("abline(h=0, lty=2)")
	}
	doItAndPrint("par(mfrow=.mfrow)")
	logger("remove(.dfbetas, .mfrow)")
	remove(.dfbetas, .mfrow, envir=.GlobalEnv)
}

CoxDfbeta <- function(){ # works for survreg models as well
	command <- paste(".dfbeta <- as.matrix(residuals(", ActiveModel(), ', type="dfbeta"))', sep="")
	doItAndPrint(command)
	if (coxphP()){
		command <- paste("colnames(.dfbeta) <- names(coef(", ActiveModel(), "))", sep="")
		doItAndPrint(command)
	}
	ncol <- ncol(.dfbeta)
	doItAndPrint(paste(".mfrow <- par(mfrow = mfrow(", ncol, "))", sep=""))
	for (col in colnames(.dfbeta)){
		doItAndPrint(paste('plot(.dfbeta[,"', col, '"])', sep=""))
		doItAndPrint("abline(h=0, lty=2)")
	}
	doItAndPrint("par(mfrow=.mfrow)")
	logger("remove(.dfbeta, .mfrow)")
	remove(.dfbeta, .mfrow, envir=.GlobalEnv)
}

MartingalePlots <- function(){
	.activeModel <- ActiveModel()
	command <- paste(".NullModel <- update(", .activeModel, ", . ~ 1)", sep="")
	doItAndPrint(command)
	doItAndPrint('.residuals <- residuals(.NullModel, type="martingale")')
	command <- paste(".X <- padNA(model.matrix(", .activeModel, 
		"), residuals(", .activeModel, "))", sep="")
	doItAndPrint(command)
	coefs <- names(coef(eval(parse(text=.activeModel))))
	ncoef <- length(coefs)
	doItAndPrint(paste(".mfrow <- par(mfrow = mfrow(", ncoef, "))", sep=""))
	activeDataSet <- ActiveDataSet()
	for (coef in coefs){
		x <- paste('.X[,"', coef, '"]', sep="")
		command <- if (length(unique(eval(parse(text=x)))) < 10)
				paste("plot(", x, ', .residuals, xlab="', coef,
					'", ylab="Martingale residuals from null model")', sep="")
			else paste("scatter.smooth(", x, ', .residuals, xlab="', coef,
					'", ylab="Martingale residuals from null model", family="gaussian")', sep="")
		doItAndPrint(command)
		doItAndPrint(paste("abline(lm(.residuals ~ ", x, "), lty=2)", sep=""))
	}
	doItAndPrint("par(mfrow=.mfrow)")
	logger("remove(.NullModel, .residuals, .X, .mfrow)")
	remove(.residuals, .X, .mfrow, envir=.GlobalEnv)	
}

PartialResPlots <- function(){
	command <- paste(".residuals <- residuals(", ActiveModel(), ', type="partial")', sep="")
	doItAndPrint(command)
	command <- paste(".fitted <- predict(", ActiveModel(), ', type="terms")', sep="")
	doItAndPrint(command)
	terms <- colnames(.residuals)
	nterms <- length(terms)
	doItAndPrint(paste(".mfrow <- par(mfrow = mfrow(", nterms, "))", sep=""))
	activeDataSet <- ActiveDataSet()
	for (term in terms){
		x <- paste('.fitted[,"', term, '"]', sep="")
		command <- if (length(unique(eval(parse(text=x)))) < 10)
				paste("plot(", x, ', .residuals[,"', term, '"], xlab="', term,
					'", ylab="Partial residuals")', sep="")
			else paste("scatter.smooth(", x, ', .residuals[,"', term, '"], xlab="',
					term, '", ylab="Partial residuals", family="gaussian")', sep="")
		doItAndPrint(command)
		command <- paste("abline(lm(", '.residuals[,"', term, '"] ~ ', x, "))", sep="")
		doItAndPrint(command)
	}
	doItAndPrint("par(mfrow=.mfrow)")
	logger("remove(.residuals, .fitted, .mfrow)")
	remove(.residuals, .fitted, .mfrow, envir=.GlobalEnv)	
}
