cv.AWEnet <-
function (X, Y, delta, weight, lambda2, maxit, K = 10, fraction = seq(from = 0, to = 1, length = 100), 
            plot.it = F, se = TRUE, AEnet=T, all.folds=NULL) 
{
  #require(lars)
   bds<-sort(lambda2)
	cv.Enet <-
	function (X, Y, delta, weight, lambda2, maxit, K, fraction, AEnet) 
	{
  	if (is.null(all.folds))
    	all.folds <- cv.folds(length(Y), K)
  	residmat <- matrix(0, length(fraction), K)
  	for (i in seq(K)) {
    	omit <- all.folds[[i]]
    	if (AEnet)
      fit <- AEnet.aft(X, Y, delta, weight, lambda2, maxit)
    	else
      fit <- WEnet.aft(X, Y, delta, weight, lambda2, maxit)
    	fit <- predict(fit, X[omit,, drop = FALSE], mode = "fraction",s = fraction)$fit
    	if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    	residmat[, i] <- apply((Y[omit] - fit)^2, 2, mean)
  	}
  	cv <- apply(residmat, 1, mean)
  	cv.error <- sqrt(apply(residmat, 1, var)/K)
  	object <- list(index = fraction, cv = cv, cv.error =  cv.error,all.folds=all.folds, mode="fraction")
  	if (plot.it) 
    	plotCVLars(object, se = se)
  	invisible(object)
		}

   index<-NULL
	for (lambda in bds)
	{
      if (AEnet){
	cvEnet<-cv.Enet(X, Y, delta, weight, lambda, maxit, K , fraction, AEnet=T) 
	s<-cvEnet$index[which.min(cvEnet$cv)]
	cv.mse<-which.min(cvEnet$cv)
	cv.error<-which.min(cvEnet$cv.error)
	index<-rbind(index, c(lambda, s, cv.mse ,cv.error))}
	else{
	cvEnet<-cv.Enet(X, Y, delta, weight, lambda, maxit, K , fraction, AEnet=F) 
	gama<-cvEnet$index[which.min(cvEnet$cv)]
	s<-gama*sqrt(1+lambda)
	cv.mse<-which.min(cvEnet$cv)
	cv.error<-which.min(cvEnet$cv.error)
	index<-rbind(index, c(lambda, s, cv.mse ,cv.error))}
	}
list(index=index)
}
