predict.Bayesthresh <- function(object, ...)
{
				if(!inherits(object, "Bayesthresh"))
								stop("Use an object of class Bayesthresh")
				fixed <- coef(object)[,1]
				comp <- sum(object$compVar[,1])
				if(object$link == 'Gaussian'){
								if(is.null(object$X)){
												x <- rep(1,length(object$Y))
												if(length(object$fl)==1){
																rad <- random.effects(object)[,1]
												}
												pred <- x*fixed+object$Zl%*%rad
								}
								else
								{
												x <- object$X
												ifelse(length(object$compVar[,1])==2,rad <- random.effects(object)[,1],
																			rad <- unlist(lapply(random.effects(object), function(x)x[,1])))
								}
								thres <- c(0, object$Ntau,Inf)
								pred <- x%*%fixed+object$Zl%*%rad
								ty <- rep(0,length(object$Y))
								for(i in 2:(length(thres)-1)) {
												ty <- ty+(pnorm(thres[i+1]-pred,0,sqrt(comp))-pnorm(thres[i]-pred,0,sqrt(comp)))*(i+1)
								}
				}
				else{
								if(object$link != 't') stop('Link function is not recognized')
								if(is.null(object$X)){
												x <- rep(1,length(object$Y))
												if(length(object$fl)==1){
																rad <- random.effects(object)[,1]
												}
												pred <- x*fixed+object$Zl%*%rad
								}
								else
								{
												x <- object$X
												ifelse(length(object$compVar[,1])==2,rad <- random.effects(object)[,1],
																			rad <- unlist(lapply(random.effects(object), function(x)x[,1])))
								}
								thres <- c(0, object$Ntau,Inf)
								dfg <- length(object$Y)-rankMatrix(cbind(object$X,object$W))
								pred <- x%*%fixed+object$Zl%*%rad
								ty <- rep(0,length(object$Y))
								for(i in 2:(length(thres)-1)) {
												ty <- ty+(pt(thres[i+1]-pred,dfg)-pt(thres[i]-pred,dfg))*(i+1)
								}
				}
				return(ty)
}

