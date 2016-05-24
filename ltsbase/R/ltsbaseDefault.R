ltsbaseDefault <-
function(xdata,y,alpha=alpha, by=by){

	
	est=ltsbase(xdata,y,alpha=alpha, by=by)
        est$fitted.values <- as.matrix(xdata)%*%as.matrix(est$list.coef.all)

	est$residuals <- y - est$fitted.values
        est$call=match.call()

class(est)= "ltsbase"

results=list(fitted.val=est$fitted.values,res=est$residuals)
results
}
