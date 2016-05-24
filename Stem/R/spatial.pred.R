`spatial.pred` <-
function(mu1,mu2,Sigma11,Sigma12,Sigma22,X1){
	pred = mu2 + Sigma12 %*% solve(Sigma11) %*% (X1-mu1)
	var.pred = Sigma22 - Sigma12 %*%solve(Sigma11)%*%t(Sigma12)
	se.pred = sqrt(diag(var.pred))
	return(list(pred=pred,se.pred=se.pred))
}

