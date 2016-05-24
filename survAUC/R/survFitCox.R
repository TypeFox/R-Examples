
################################################################
###			Survivor-Function of Cox-Modell					####
################################################################
## lp			- the vector of linear predictor of training data
## response		- the Surv(.,.) Outcome of training data
## lpnew		- the vector of linear predictors of test data

survFit.COX <- function(lp, response, lpnew)
{
	time <- response[,1]
	event <- response[,2]
	n_time <- length(time)
	n_lp <- length(lp)
	n_lpnew <- length(lpnew)
	
	ans <- .Call("survfit_cox", 
				 as.numeric(lp), 
				 as.numeric(time), 
				 as.numeric(event), 
				 as.integer(n_time),
				 as.integer(n_lp),
				 as.numeric(lpnew),
				 as.integer(n_lpnew),
				 PACKAGE="survAUC")
	ans
}




################################################################
###						Partial Likelihood					 ###
###						  Cox-Model  						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors of training data

PartLCox <- function(Surv.rsp, lp){
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	n.event <- length(event)
	n.lp <- length(lp)
	if(n.lp == 1)
	lp <- rep(lp,n.event)
	
	ans <- .C("partLCox",
			  as.numeric(stime),
			  as.numeric(event),
			  as.integer(n.event),
			  as.numeric(lp),
			  as.integer(n.lp),
			  as.numeric(0),
			  PACKAGE="survAUC")
	ans[[6]]
}



################################################################
###						censoring weights					 ###
###						  COX-Model 						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## Surv.rsp.new	- Surv(.,.) Outcome of test data
## times		- time points

censWeights <- function(Surv.rsp, Surv.rsp.new, times)
{
## Surv-train
	stime <- Surv.rsp[,1]
	event <- 1-Surv.rsp[,2]
	
## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]
	
	n.stime <- length(stime)
	n.stime.new <- length(stime.new)
	n.times <- length(times)
	weights <- matrix(0, ncol=n.stime.new, nrow=n.times)
	
	ans <- .C("cens_weights",
			  as.numeric(times),
			  as.integer(n.times),
			  as.numeric(stime),
			  as.numeric(event),
			  as.integer(n.stime),
			  as.numeric(stime.new),
			  as.numeric(event.new),
			  as.integer(n.stime.new),
			  as.numeric(weights),
			  PACKAGE="survAUC")
	ans
}







################################################################
###						Partial Likelihood					 ###
###						  INDV								 ###
################################################################
## stime		- Surv(.,.)[,1] Outcome of survival response
## time			- time point
## lp			- vector of linear predictors of training data

PartLCoxIndiv <- function(stime, time, lp){
	
	ans <- .C("partLCoxIndiv",
			  as.numeric(stime),
			  as.numeric(time),
			  as.integer(length(stime)),
			  as.numeric(lp),
			  as.numeric(rep(0,length(stime))),
			  PACKAGE="survAUC")
	ans
}


