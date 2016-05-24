theta.ini <-
function(z, X, CorModels, use.nugget,
	use.anisotropy, dist.hydro.data, x.dat, y.dat, REs)
{
	n.models <- length(CorModels)
	var.resid <- mean((z - X %*% mginv(t(X) %*% X) %*% t(X) %*% z)^2)
	theta <- NULL
    scale <- NULL ## What scale is the parameter on - to simplify back transformation
    type <- NULL
    terms <- NULL
	if(length(grep("tailup",CorModels)) > 0){
		if(length(grep("tailup",CorModels)) > 1)
			stop("Cannot have more than 1 tailup model")
		theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
			log(mean(dist.hydro.data))),ncol = 1))
		scale <- c(scale,c("log","log"))
        type <- c(type,c("parsill","range"))
        terms <- c(terms,rep(CorModels[grep("tailup",CorModels)],
			times = 2))
	}
	if(length(grep("taildown",CorModels)) > 0){
		if(length(grep("taildown",CorModels)) > 1)
			stop("Cannot have more than 1 taildown model")
		theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
			log(mean(dist.hydro.data))),ncol = 1))
		scale <- c(scale,c("log","log"))
        type <- c(type,c("parsill","range"))
        terms <- c(terms,rep(CorModels[grep("taildown",CorModels)],
			times = 2))
	}
	if(length(grep("Euclid",CorModels)) > 0){
		if(length(grep("Euclid",CorModels)) > 1)
			stop("Cannot have more than 1 Euclidean model")
		dist.Euclid.data <- distGeo(x.dat,y.dat,x.dat,y.dat)
		if(use.anisotropy == FALSE) {
        	theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
				log(mean(dist.Euclid.data))),ncol = 1))
			scale <- c(scale,c("log","log"))
        	type <- c(type,c("parsill","range"))
        	terms <- c(terms,rep(CorModels[grep("Euclid",CorModels)],
				times = 2))
		}
		else {
			theta <- rbind(theta,matrix(c(log(.9/n.models*var.resid),
				log(mean(dist.Euclid.data)),0,0),ncol = 1))
			scale <- c(scale,c("log","log","logistic","logistic180"))
        	type <- c(type,c("parsill","range","axratio","rotate"))
        	terms <- c(terms,rep(CorModels[grep("Euclid",CorModels)],
				times = 4))
		}
	}
	if(length(REs) > 0) {
            theta <- rbind(theta,matrix(rep(log(.9/n.models*var.resid),
                                            times = length(REs)), ncol = 1))
            scale <- c(scale,rep("log", times = length(REs)))
            type <- c(type,rep("parsill", times = length(REs)))
            terms <- c(terms, names(REs))
	}
	if(use.nugget == TRUE) {
        if(is.null(theta)) {
            theta <- log(var.resid)
            scale <- "log"
            type <- "parsill"
            terms <- "Nugget"
        } else {
            theta <- c(theta,log(.1*var.resid))
            scale <- c(scale,"log")
            type <- c(type,"parsill")
            terms <- c(terms,"Nugget")
        }
	}

    attr(theta,"scale") <- scale
    attr(theta,"type") <- type
    attr(theta,"terms") <- terms

	theta

}

