### Return the sample of MCMC process

MCMCsample <- function(object){
				if(!inherits(object, "Bayesthresh"))
								stop("Use an object of class Bayesthresh")
				if(object$Write!=TRUE) stop("'Write' is not TRUE")
				theta <- data.frame(object$outp.mcmc[[1]])
				colnames(theta) <- object$NamesTheta
				cutpoints <- object$outp.mcmc[[2]]
				Variance <- data.frame(object$outp.mcmc[[3]])
				colnames(Variance) <- object$NomesVu
				list(Theta=theta, Variance=Variance, cutpoints=cutpoints)

}
