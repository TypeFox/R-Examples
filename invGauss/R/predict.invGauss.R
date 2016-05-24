predict.invGauss <- function(object, newdata, warn = T, ...){
##
## PREDICTS SURVIVAL FOR EACH LINE IN A NEW DATA FILE

.coef <- object$coef
.link.mu <- object$link.mu
.link.c <- object$link.c

.tau <- .coef["tau"]
.intercept.mu <- .coef["Intercept.mu"] # NOTE: ASSUMES INTERCEPT IS ALWAYS IN MODEL
.intercept.c <- .coef["Intercept.c"] # NOTE: ASSUMES INTERCEPT IS ALWAYS IN MODEL

.betas <- .coef[-match(c("tau", "Intercept.mu", "Intercept.c"), names(.coef))]

.mu.c <- substring(names(.betas), first = nchar(names(.betas)))


.betas.mu <- .betas[.mu.c == "u"]
.betas.c <- .betas[.mu.c == "c"]

.match.mu <- match(names(.betas.mu), paste(names(newdata), ".mu", sep = "")) # FIND POSITIONS OF RELEVANT COVARIATES
.match.c <- match(names(.betas.c), paste(names(newdata), ".c", sep = "")) # FIND POSITIONS OF RELEVANT COVARIATES


.X.mu <- newdata[, .match.mu, drop = F] # EXTRACT DATA WITH ONLY RELEVANT COVARIATES
.X.c <- newdata[, .match.c, drop = F] # EXTRACT DATA WITH ONLY RELEVANT COVARIATES


.mu <- .link.mu(as.matrix(.X.mu) %*% .betas.mu + .intercept.mu) # COMPUTES EACH MU VALUE
.c1 <- .link.c(as.matrix(.X.c) %*% .betas.c + .intercept.c) # COMPUTES EACH STARTING VALUE



if(warn) cat("BE CAREFUL: the prediction function is a beta version, it requires the time variable to be names, just that, 'time'. It also assumes intercepts are present. Use at own risk...\n") 

.input <- data.frame(t = newdata[,"time"], tau = .tau, mu = .mu, c1 = .c1, row.names = NULL)


.fb <- apply(.input, 1, function(x) f.fb(x[1], x[3], x[2], x[4])) # WATCH ARGUMENT SEQUENCE!
.B <- apply(.input, 1, function(x) f.B(x[1], x[3], x[2], x[4])) # WATCH ARGUMENT SEQUENCE!

.ut <- data.frame(newdata, mu = .mu, c = .c1, dens = .fb, surv = .B)

return(.ut)

return(.betas)


.datanames <- names(newdata)






if(length(.coef) == 3){
	cat("\nNote: no covariates used in model\n\n")
}

else{







	
}



if(F){
# EXTRACT time, status AND covariates
	model <- model.frame(formula = object$call$formula, data = newdata)
	response <- model.extract(model, "response")
	time <- response[, "time"]
###	status <- response[, "status"]
###	assign("no.censor", all(status == 1), f = 1)	#
###	assign("no.censor", all(status == 1), env = .GlobalEnv)	#
###	if(all(status != 1))
###		stop("No events!")	#
	covar <- model.matrix(formula, data = data)
###	assign("ncovar", dim(covar)[2], f = 1)
###	assign("ncovar", dim(covar)[2], env = .GlobalEnv)
	termnames <- c("mu", "tau", dimnames(covar)[[2]])	#
}


### return(cbind(.X, .input))

}
