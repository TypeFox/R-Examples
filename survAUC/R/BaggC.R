

##############################################################
## C-Statistic suggest by Begg
##############################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data


BeggC <- function(Surv.rsp, Surv.rsp.new, lp, lpnew){

## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]
## Times
	times <- stime.new
	n_times <- length(times)
	n_stime <- length(stime)
	n_lp <- length(lp)
	n_stime_new <- length(stime.new)
	n_lpnew <- length(lpnew)
	if(n_stime != n_lp)
		stop(" 'Surv.rsp' and 'linear predictors' must have the same length!\n")
	if(n_stime_new != n_lpnew)
		stop(" 'Surv.rsp.new' and 'linear predictors new' must have the same length!\n")

	#### Cox survival function estimates for lpnew	
	surv.cox <- .Call("survfit_cox",
					  as.numeric(lp), 
					  as.numeric(stime),
					  as.numeric(event), 
					  as.integer(n_stime),
					  as.integer(n_lp),
					  as.numeric(lpnew),
					  as.integer(n_lpnew),
					  PACKAGE="survAUC")
	#### C-Statistic
	c.begg <- .C("c_begg",
				   as.numeric(stime.new), 
				   as.numeric(event.new),
				   as.integer(n_stime_new),
				   as.numeric(times),
				   as.integer(n_times),
				   as.numeric(lp),
				   as.numeric(lpnew),
				   as.numeric(surv.cox[[1]]),
				   as.numeric(surv.cox[[2]]),
				   as.integer(length(surv.cox[[2]])),
				   as.numeric(vector("numeric",length=1)),
				   PACKAGE="survAUC")
	c.begg[[11]]
}


