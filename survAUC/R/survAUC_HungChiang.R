
################################################################
###						Hung & Chiang						####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times


AUC.hc <- function(Surv.rsp, Surv.rsp.new, lpnew, times)
{
## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]

	n_time <- length(times)
	n_stime <- length(stime)
	n_stime_new <- length(stime.new)
	n_lpnew <- length(lpnew)
	auc <- vector("numeric",n_time)
	
	ans <- .C("Hung_Chiang",
			  as.numeric(times),
			  as.integer(n_time),
			  as.numeric(stime),
			  as.numeric(event),
			  as.integer(n_stime),
			  as.numeric(stime.new),
			  as.numeric(event.new),
			  as.integer(n_stime_new),
			  as.numeric(lpnew),
			  as.integer(n_lpnew),
			  as.numeric(auc),
			  as.numeric(0),
			  PACKAGE="survAUC")
	erg <- list(auc=ans[[11]], times=times, iauc=ans[[12]])
	class(erg) <- "survAUC"
	erg
}


