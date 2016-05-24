

################################################################
###						Song & Zhou			     			####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times


sens.sh <- function(Surv.rsp, lp, lpnew, times, type="incident"){

	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
	n_stime <- length(stime)
	n_lp <- length(lp)
	n_lpnew <- length(lpnew)
	
	type_sens <- charmatch(type,c("incident","cumulative"))
	if(is.na(type_sens))
		stop("\nThe value of 'type' must be one of 'incident' or 'cumulative'\n")
	
	ans <- .Call("sens_SZ",
				 as.numeric(unique(sort(lpnew))),
				 as.numeric(times),
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n_stime),
				 as.numeric(lp),
				 as.integer(n_lp),
				 as.numeric(lpnew),
				 as.integer(n_lpnew),
				 as.logical(type_sens-1),
				 PACKAGE="survAUC")
	ans
}






################################################################
###						Song & Zhou						####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times


spec.sh <- function(Surv.rsp, lp, lpnew, times){

	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
	n_stime <- length(stime)
	n_lp <- length(lp)
	n_lpnew <- length(lpnew)
	
	ans <- .Call("spez_SZ",
				 as.numeric(unique(sort(lpnew))),
				 as.numeric(times),
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n_stime),
				 as.numeric(lp),
				 as.integer(n_lp),
				 as.numeric(lpnew),
				 as.integer(n_lpnew),
				 PACKAGE="survAUC")
	ans
}




################################################################
###						Song & Zhou						####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times


AUC.sh <- function(Surv.rsp, Surv.rsp.new=NULL, lp, lpnew, times, type="incident", savesensspec=FALSE)
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
	type_sens <- charmatch(type,c("incident","cumulative"))
	if(is.na(type_sens))
		stop("\nThe value of 'type' must be one of 'incident' or 'cumulative'\n")
	n_stime <- length(stime)
	n_stime.new <- length(stime.new)
	n_lp <- length(lp)
	n_lpnew <- length(lpnew)
	
	ans <- .Call("auc_SZ",
				 as.numeric(unique(sort(lpnew))),
				 as.numeric(times),
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n_stime),
				 as.numeric(stime.new),
				 as.numeric(event.new),
				 as.integer(n_stime.new),
				 as.numeric(lp),
				 as.integer(n_lp),
				 as.numeric(lpnew),
				 as.integer(n_lpnew),
				 as.logical(type_sens-1),
				 PACKAGE="survAUC")
	if(!savesensspec){
		erg <- list(auc=ans[[1]], times=ans[[2]], iauc=ans[[5]])
	}else{
		erg <- ans
	}
	class(erg) <- "survAUC"
	erg
}
