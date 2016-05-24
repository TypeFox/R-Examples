testCrossValidation <- function(model,Kfold=c(2,5,10,20,30,40,dim(model$data$X)[1]),N=10){
	
	n <- dim(model$data$X)[1]

	out1 <- matrix(0,ncol=N,nrow=length(Kfold))	
	out2 <- matrix(0,ncol=N,nrow=length(Kfold))	

	for(i in 1:length(Kfold)){
	  if(Kfold[i]==n){
	    tmp_ <- crossValidation(model,K=Kfold[i])
		out1[i,1:N] <- tmp_$Q2
		out2[i,1:N] <- tmp_$RMSE_CV
	  } else for (j in 1:N){
		tmp_ <- crossValidation(model,K=Kfold[i])
		out1[i,j] <- tmp_$Q2
		out2[i,j] <- tmp_$RMSE_CV
	  }
	}

	plot(rep(Kfold,each=N),as.numeric(t(out1)),type='p',xlab='Number of folds k',ylab='Q2',main='Q2 criterion obtained by k-fold cross validation')
	par(ask=TRUE)
	plot(rep(Kfold,each=N),as.numeric(t(out2)),type='p',xlab='Number of folds k',ylab='RMSE_CV',main='RMSE criterion obtained by k-fold cross validation')
	par(ask=FALSE)

	return(list(r2=out1,rmse=out2))
}