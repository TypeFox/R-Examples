###############################################################################
#' SMS-EMOA: S-Metric-Selection Evolutionary Multi-objective Optimization Algorithm
#' 
#' Straight forward SMS-EMOA implementation. This function is used to optimize several surrogate models
#' when doing multi objective optimization with SPOT. See: \code{\link{spotParetoOptMulti}}.
#'
#' @param f target function to be optimized of  type f(x)=y where both x and y are vectors.
#' 			The target function should return a vector of length(y) containing NAs if the input vector x contains NA values.
#' @param lower the lower boundary vector of the decision space
#' @param upper the upper boundary vector of the decision space
#' @param ... further settings relayed to \code{f}
#' @param control list of parameters (defaults are: mu=100L, sbx.n=15, sbx.p=0.7, pm.n=25, pm.p=0.3)
#' @return list with archive of solutions, active Pareto front and others
#' @author O. Mersmann
#'
#' @references N. Beume, B. Naujoks, and M.Emmerich. \emph{SMS-EMOA: Multi-objective selection based on dominated hypervolume}.
#' European Journal of Operational Research, 181(3):1653--1669, 2007. \cr \cr
#' Link to the sms_emoa code by Olaf Mersmann: \url{http://git.p-value.net/p/emoa.git/tree/examples/sms_emoa.r} 
#' 
#' @export
###############################################################################
spotSmsEmoa <- function(f, lower, upper, ...,
                     control=list(mu=100L,
                       sbx.n=15, sbx.p=0.7,
                       pm.n=25, pm.p=0.3
                       )) {
  ## Extract control parameters:
  default <- formals(sys.function())$control
  control <- steady_state_emoa_control(f, lower, upper, ..., control=control, default=default)
  control <- sbx_control(f, upper, lower, ..., control=control, default=default)
  control <- pm_control(f, upper, lower, ..., control=control, default=default)  
  control$ref <- emoa::coalesce(control[["ref"]], rep(11, control$d))

  ## Tracking variables:
  X <- matrix(0, nrow=control$n, ncol=control$maxeval)
  Y <- matrix(0, nrow=control$d, ncol=control$maxeval)
  dob <- rep(-1L, control$maxeval)
  eol <- rep(-1L, control$maxeval)
  
  ## Random inital population:
  X[, 1:control$mu] <- replicate(control$mu, runif(control$n, lower, upper))
  Y[, 1:control$mu] <- sapply(1:control$mu, function(i) f(X[,i],...))

  neval <- control$mu       ## Count the number of function evaluations
  active <- 1:control$mu    ## Indices of individuals that are in the current pop.

  ## Save some common control parameters into the current
  ## environment. This saves a few msec of execution time...
  crossover <- control$crossover
  mutate <- control$mutate
  maxeval <- control$maxeval
  #logger <- control$logger
  
  #logger$start("sms_emoa")
  while(neval < maxeval) {
    ############################################################
    ## Variation:
    parents <- sample(active, 2)
    child <- crossover(X[, parents])[,sample(c(1, 2), 1)]
    x <- mutate(child)

    ## Add new individual:
    neval <- neval + 1
    X[, neval] <- x
    Y[, neval] <- f(x,...)
    dob[neval] <- neval ## For a steady state emoa this is trivial...
    active <- c(active, neval)

    ############################################################
    ## Selection:
    i <- nds_hv_selection(Y[, active])

    ## Remove the i-th active individual:
    eol[active[i]] <- neval
    active <- active[-i]
  }  
  structure(list(X=X, Y=Y,
                        dob=dob,
                        eol=eol,
                        par=X[,active], value=Y[,active]),
                   class="emoa_result")
}


###############################################################################
#' SMS-EMOA: S-Metric-Selection Evolutionary Multi-objective Optimization Algorithm
#' 
#' Straight forward SMS-EMOA implementation, but supported by a Kriging model when selecting new individuals.
#'
#' @param f target function to be optimized of  type f(x)=y where both x and y are vectors.
#' 			The target function should return a vector of length(y) containing NAs if the input vector x contains NA values.
#' @param lower the lower boundary vector of the decision space
#' @param upper the upper boundary vector of the decision space
#' @param ... further settings relayed to \code{f}
#' @param control list of parameters (defaults are: mu=100L, sbx.n=15, sbx.p=0.7, pm.n=25, pm.p=0.3)
#' @return list with archive of solutions, active Pareto front and others
#'
#' @seealso N. Beume, B. Naujoks, and M.Emmerich. \emph{SMS-EMOA: Multiobjective selection based on dominated hypervolume}.
#' European Journal of Operational Research, 181(3):1653--1669, 2007. \cr \cr
#' Link to the sms_emoa code by Olaf Mersmann: \url{http://git.p-value.net/p/emoa.git/tree/examples/sms_emoa.r} 
#' 
#' @export
###############################################################################
spotSmsEmoaKriging <- function(f, lower, upper, ...,
                     control=list(mu=100L,
                       sbx.n=15, sbx.p=0.7,
                       pm.n=25, pm.p=0.3
                       )) {
  ## Extract control parameters:
  default <- formals(sys.function())$control
  control <- steady_state_emoa_control(f, lower, upper, ..., control=control, default=default)
  control <- sbx_control(f, upper, lower, ..., control=control, default=default)
  control <- pm_control(f, upper, lower, ..., control=control, default=default)  
  control$ref <- emoa::coalesce(control[["ref"]], rep(11, control$d))

  ## Tracking variables:
  X <- matrix(0, nrow=control$n, ncol=control$maxeval)
  Y <- matrix(0, nrow=control$d, ncol=control$maxeval)
  dob <- rep(-1L, control$maxeval)
  eol <- rep(-1L, control$maxeval)
  
  ## Random inital population:
  X[, 1:control$mu] <- replicate(control$mu, runif(control$n, lower, upper))
  Y[, 1:control$mu] <- sapply(1:control$mu, function(i) f(X[,i],...))

  neval <- control$mu       ## Count the number of function evaluations
  active <- 1:control$mu    ## Indices of individuals that are in the current pop.

  ## Save some common control parameters into the current
  ## environment. This saves a few msec of execution time...
  crossover <- control$crossover
  mutate <- control$mutate
  maxeval <- control$maxeval
  #logger <- control$logger
  
  #logger$start("sms_emoa")
  while(neval < maxeval) {
		#build minimal config to use spot kriging predictor
		rawB=data.frame(cbind(t(X[,1:neval]),t(Y[,1:neval])))
		spotConfig<-list(alg.roi= matrix(0,nrow=nrow(X),ncol=3),
			io.verbosity=0,
			seq.infill=NA,
			mco.refPoint=apply(Y,1,max)+1,
			alg.resultColumn = paste(rep("y",nrow(Y)),1:nrow(Y),sep="")
		)
		rownames(spotConfig$alg.roi)=paste(rep("x",nrow(X)),1:nrow(X),sep="")
		colnames(rawB)=c(rownames(spotConfig$alg.roi),spotConfig$alg.resultColumn)
		spotConfig<-spotPredictForrester(rawB,NULL,t(X[,1:2]),spotConfig)

		getFitness <- function(x){
		if(any(is.na(x))){rep(NA,length(spotConfig$alg.resultColumn))}  #sms-emoa starts with NA values to determine dimension
		as.numeric(eval(call("spotPredictForrester"
						, NULL 
						, NULL
						, as.data.frame(x)
						, spotConfig
						, spotConfig$seq.modelFit #external fit is used, model is only evaluated not build, therefore the NULLS are no prob
						))$seq.largeDesignY);				
		}

		nn=control$npoints
		x=matrix(0,nrow=nrow(X),ncol=nn)
		y=NULL
		for(i in 1:nn){
		parents <- sample(active, 2)
		child <- crossover(X[, parents])[,sample(c(1, 2), 1)]
		x[,i] <- mutate(child)
		y <- rbind(y,dominated_hypervolume(cbind(Y[,active],getFitness(x[,i])))-dominated_hypervolume(Y[,active]))		
		}
		x<-x[,which.max(y)]
		## Add new individual:
		neval <- neval + 1
		X[, neval] <- x
		Y[, neval] <- f(x,...)
		dob[neval] <- neval ## For a steady state emoa this is trivial...
		active <- c(active, neval)
		############################################################
		## Selection:
		i <- nds_hv_selection(Y[, active])

		## Remove the i-th active individual:
		eol[active[i]] <- neval
		active <- active[-i]
  }  
	structure(list(X=X, Y=Y,
                        dob=dob,
                        eol=eol,
                        par=X[,active], value=Y[,active]),
                   class="emoa_result")
 }

  