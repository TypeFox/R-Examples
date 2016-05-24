oraclefit <-
function(data,testdata,family="gaussian"){
	fit = NULL
	x = data$x0
	y = data$y
	xtest = testdata$x0
	ytest = testdata$y
	#x = scale(x)
	
	lmfit = glm(y~x,family=family)
	fit$a0 = coef(lmfit)[1]
	fit$beta = coef(lmfit)[-1]
	#fit$beta = fit$beta/attr(x,'scaled:scale')
	#fit$a0 = fit$a0-sum(fit$beta*attr(x,'scaled:center'))
	
	fit$mse = sum((fit$beta-data$beta0)^2)+fit$a0^2
	yfit = cbind(1,xtest)%*%c(fit$a0,fit$beta)
	if(family=='gaussian')
	fit$Rsq = 1-sum((ytest-yfit)^2)/sum((ytest-mean(ytest))^2)
	else fit$error = mean(ytest!=(yfit>0))
	return(fit=fit)
}
