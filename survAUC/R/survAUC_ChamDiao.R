
################################################################
###						Chambers & Diao						####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times


AUC.cd <- function(Surv.rsp, Surv.rsp.new = NULL, lp, lpnew, times)
{
## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
## Surv-test
	if(!is.null(Surv.rsp.new)){
		stime.new <- Surv.rsp.new[,1]
		event.new <- Surv.rsp.new[,2]
	}else{
		stime.new <- NULL
		event.new <- NULL
	}
	n_stime <- length(stime)
	n_lp <- length(lp)
	n_lpnew <- length(lpnew)
	
	erg <- .Call("Cham_Diao",
				 as.numeric(lp),
				 as.numeric(times),
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n_stime),
				 as.numeric(stime.new),
				 as.numeric(event.new),
				 as.integer(length(stime.new)),
				 as.integer(n_lp),
				 as.numeric(lpnew),
				 as.integer(n_lpnew),
				 PACKAGE="survAUC")
	class(erg) <- "survAUC"
	erg
}


