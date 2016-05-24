##############################################################
## Uno sensetifity 
##############################################################
# Surv.rsp = Zielvariable train, Surv-Objekt (time, status)
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praedikt. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll

sens.uno <- function(Surv.rsp, Surv.rsp.new, lpnew, times){
	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	ERG <- .C("sens_uno",
			  as.numeric(rep(1, n_t*(n_th+1))),
			  as.numeric(Surv.rsp[,1]),
			  as.numeric((1-Surv.rsp[,2])),
			  as.numeric(thresh),
			  as.numeric(times),
			  as.numeric(lpnew),
			  as.numeric(Surv.rsp.new[,1]),
			  as.numeric(Surv.rsp.new[,2]),
			  as.integer(n_th),
			  as.integer(n_t),
			  as.integer(dim(Surv.rsp.new)[1]),
			  as.integer(dim(Surv.rsp)[1]),
			  PACKAGE="survAUC")
	matrix(ERG[[1]], n_t, n_th+1)
}


##############################################################
## Uno specificity
##############################################################
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praedikt. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll

spec.uno <- function(Surv.rsp.new, lpnew, times){
	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	ERG <- .C("spec_uno", 
			  as.numeric(rep(0, n_t*(n_th+1))), 
			  as.numeric(thresh), 
			  as.numeric(times),
			  as.numeric(lpnew), 
			  as.numeric(Surv.rsp.new[,1]), 
			  as.integer(n_th),
			  as.integer(n_t), 
			  as.integer(dim(Surv.rsp.new)[1]),
			  PACKAGE="survAUC")
	matrix(ERG[[1]], n_t, n_th+1)
}



##############################################################
## Uno AUC
##############################################################
# Surv.rsp = Zielvariable train, Surv-Objekt (time, status)
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praediktoren aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll
# weight = Welche Gewichtung der Integrated AUC?; rescale oder conditional.



AUC.uno <- function(Surv.rsp, Surv.rsp.new, lpnew, times, savesensspec=FALSE){

	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	
	#### Sensetivity, Specificity and AUC.
	auc.uno <- .C("auc_uno",
				  as.numeric(vector("numeric",length=n_t)),
				  as.numeric(0),
				  as.numeric(vector("numeric",length=n_t*(n_th+1))+1),
				  as.numeric(vector("numeric",length=n_t*(n_th+1))),
				  as.numeric(Surv.rsp[,1]),
				  as.numeric(1-Surv.rsp[,2]),
				  as.numeric(thresh), 
				  as.numeric(times),
				  as.numeric(lpnew), 
				  as.numeric(Surv.rsp.new[,1]),
				  as.numeric(Surv.rsp.new[,2]),
				  as.integer(n_th),
				  as.integer(n_t), 
				  as.integer(dim(Surv.rsp.new)[1]),
				  as.integer(dim(Surv.rsp)[1]),
				  PACKAGE="survAUC")
	if(!savesensspec){
		erg <- list(auc=auc.uno[[1]], times=auc.uno[[8]], iauc=auc.uno[[2]])
	}else{
		erg <- list(auc=auc.uno[[1]], times=auc.uno[[8]], iauc=auc.uno[[2]],
			 sens=matrix(auc.uno[[3]], n_t, n_th+1), 
			 spec=matrix(auc.uno[[4]], n_t, n_th+1),
			 thresh=auc.uno[[7]])
	}
	class(erg) <- "survAUC"
	erg
}


