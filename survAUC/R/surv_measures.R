


################################################################
###						Gonen and Hellers					 ###
###              Concordance Index for Cox models			 ###
################################################################
## lpnew			- the vector of linear predictors of data


GHCI <- function(lpnew){
	ans <- .C("GHCI",
			  as.numeric(lpnew),
			  as.integer(length(lpnew)),
			  as.numeric(0.0),
			  PACKAGE="survAUC")
	ans[[3]]
}






################################################################
###						Prediction Error					 ###
###						robust & brier						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## Surv.rsp.new	- Surv(.,.) Outcome of test data
## lp			- vector of linear predictors of training data
## lpnew		- vector of linear predictors of test data
## times		- vector of times
## type			- kind of prediction error curve: 'brier' or 'robust'

predErr <- function(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
					type = "brier", int.type = "unweighted")
{
	type <- charmatch( type, c("brier","robust") )
	if (is.na(type))
		stop("'type' must be one of 'brier' or 'robust'")
	int.type <- charmatch( int.type, c("weighted","unweighted") )
	if (is.na(int.type))
		stop("'int.type' must be one of 'weighted' or 'unweighted'")

	## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
	## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]

	n.times <- length(times)
	n.stime <- length(stime)
	n.stime.new <- length(stime.new)
	n.lp <- length(lp)
	n.lpnew <- length(lpnew)
	
	erg <- .Call("predError",
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n.stime),
				 as.numeric(stime.new),
				 as.numeric(event.new),
				 as.integer(n.stime.new),
				 as.numeric(times),
				 as.integer(length(times)),
				 as.numeric(lp),
				 as.integer(n.lp),
				 as.numeric(lpnew),
				 as.integer(n.lpnew),
				 as.integer(type-1),
				 as.integer(int.type-1),
				 PACKAGE="survAUC")
	class(erg) <- "survErr"
	erg
}





################################################################
###				measure by O''Quigley et al. (2005)			 ###
###						  R^2_{OXS}  						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

OXS <- function(Surv.rsp, lp, lp0)
{
	
	L <- PartLCox(Surv.rsp, lp)
	L0 <- PartLCox(Surv.rsp, lp0)	
	1 - exp( - 2 * (L-L0) / sum(Surv.rsp[,2]))
}






################################################################
###				measure by Nagelkerke						 ###
###						  R^2_{N}							 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

Nagelk <- function(Surv.rsp, lp, lp0)
{	
	L <- PartLCox(Surv.rsp, lp)
	L0 <- PartLCox(Surv.rsp, lp0)
	n <- length(lp)
	(1 - exp( - 2 * (L-L0) / n)) / (1 - exp( 2 * L0 / n))
}





################################################################
###				 measure by Xu & O''Quigley					 ###
###						  R^2_{XO}							 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

XO <- function(Surv.rsp, lp, lp0){

	time <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	n <- length(time)
	n_lp <- length(lp)
	n_lp0 <- length(lp0)
	if(n != n_lp || n_lp != n_lp0 || n != n_lp0)
		stop(" 'Surv.rsp', 'linear predictors' and 'linear predictors of null-model' must have the same length!\n")
	ans <- .C("XO",
			  as.numeric(time),
			  as.numeric(event),
			  as.integer(n), 
			  as.numeric(lp),
			  as.numeric(lp0),
			  as.numeric(0),
			  PACKAGE="survAUC")
	ans[[6]]
}




