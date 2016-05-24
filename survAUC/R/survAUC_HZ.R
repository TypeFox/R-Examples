


##############################################################
##  Kaplan-Meier
##############################################################
# Stime = follow up time
# status = status indicator (0-alive, 1-dead)


KM <- function(Stime, status){
	n_time <- length(Stime)
	km <- .C("km_Daim", 
			 vector("numeric", length(Stime)), 
			 as.numeric(Stime), 
			 as.numeric(status), 
			 as.integer(n_time),
			 PACKAGE="survAUC")
	names(km) <- c("survival","times","status","n")
	km[-4]
}



##############################################################
## weighted KM
##############################################################
# surv = AUC measurements
# times = times vector
# tmax = maximale Zeitpunkt
# weight = Welche Gewichtung der Integrated AUC?; rescale oder conditional.


weightKM <- function(Stime, status, wt=NULL, entry=NULL ){
	n_time <- length(Stime)
	if( is.null(entry) )  
		entry <- vector("numeric",n_time)
	if( is.null(wt) ) 
		wt <- vector("numeric",n_time)+1
	km <- .C("km_weight", 
			 vector("numeric", n_time), 
			 as.numeric(Stime), 
			 as.numeric(status), 
			 as.numeric(wt),
			 as.numeric(entry),
			 as.integer(n_time),
			 PACKAGE="survAUC")
	names(km) <- c("survival","times","status","weights","entry","n")
	km[-6]
}



##############################################################
## CoxWeights
##############################################################
# marker = lin. Praed. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = times vector
# status = status indicator (0 - alive, 1 - dead or event)
# thresh = Zeitpunkt, an dem ausgewertet werden soll
# entry = 


cox_weights <- function(marker, time, status, thresh, entry=NULL){
	n_time <- length(time)
	if( is.null(entry) )  
		entry <- vector("numeric",n_time)
	ans <- .C("cox_weights", 
			 as.numeric(marker), 
			 as.numeric(time), 
			 as.integer(status), 
			 as.numeric(thresh),
			 as.numeric(0),
			 as.integer(n_time),
			 PACKAGE="survAUC")
	names(ans) <- c("marker","time","status","thresh","AUC","n")
	ans[-6]
}





##############################################################
## 
##############################################################
# marker = lin. Praed. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = times vector
# status = status indicator (0 - alive, 1 - dead or event)
# thresh = Zeitpunkt, an dem ausgewertet werden soll
# entry = 





