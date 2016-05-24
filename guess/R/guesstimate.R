#' guesstimate
#' @importFrom Rsolnp solnp
#' @title Calculate item level and aggregate learning
#' @param transmatrix  transition matrix returned from \code{\link{multi_transmat}}
#' @return list with two items: parameter estimates and estimates of learning
#' @export

guesstimate <- function(transmatrix=NULL) {
		
	# Initialize results mat
	nitems	<- nrow(transmatrix)
	nparams <- ifelse(ncol(transmatrix)==4, 4, 8)
	est.opt <- matrix(ncol=nitems, nrow=nparams)
	
	# effects
	effects	<- matrix(ncol=nitems, nrow=1)

	# calculating parameter estimates
	if (nparams == 4) {
		for (i in 1:nitems) {
			est.opt[,i]	 <- tryCatch(solnp(c(.3,.1,.1,.25), guess_lik, eqfun = eqn1, eqB = c(1), LB = rep(0,4), UB = rep(1,4), data=transmatrix[i,])[[1]], error=function(e) NULL)
		}

		effects[,1:nitems] <- est.opt[2,]

	} else {
		for (i in 1:nitems) {
			est.opt[,i]	 <- tryCatch(solnp(c(.3,.1,.2,.05,.1,.1,.05,.25), guessdk_lik, eqfun = eqn1dk, eqB = c(1), LB = rep(0,8), UB = rep(1,8), data=transmatrix[i,])[[1]], error=function(e) rep(NA,8))
		}
		
		effects[,1:nitems] 	<- est.opt[2,] + est.opt[6,]
	}
	
	# Assign row names
	if (nrow(est.opt) == 8){
		row.names(est.opt) 	<- c("lgg", "lgk", "lgc", "lkk", "lcg", "lck", "lcc", "gamma")
	} else {
		row.names(est.opt) 	<- c("lgg", "lgk",  "lkk", "gamma")
	}

	res <- list(param.lca=est.opt, est.learning=effects)

	return(invisible(res))
}