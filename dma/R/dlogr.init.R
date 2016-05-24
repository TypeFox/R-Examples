dlogr.init <-
function (x,y) {
# function to get the initial estimate of beta and var(beta)
	MyData <- data.frame(x,y)
	MyModel <- glm(y~x-1,data=MyData,family=binomial(link=logit))
	rm(MyData)	
	return(list(BetaHat=coefficients(MyModel),VarBetaHat=vcov(MyModel)))
}

