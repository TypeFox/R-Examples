# S3 method for coeficients

coef.Bayesthresh <- function(object, HPDinterval=FALSE, prob=0.95,...)
{
				if(!inherits(object, "Bayesthresh"))
								stop("Use an object of class Bayesthresh")
				fixed <- object$EfFixef
				if(HPDinterval==FALSE){
								return(fixed)
				}
				else
								if(object$Write != TRUE) stop("'Write' is not TRUE")
				mcmc.interval <- HPDinterval(mcmc(object$outp.mcmc[[1]][,c(1:dim(fixed)[1])]), prob=prob)
				rownames(mcmc.interval) <- rownames(fixed)
				list(Coefficients=fixed, HPD.interval=mcmc.interval)
}


