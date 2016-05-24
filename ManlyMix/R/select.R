
forward <- function(X, reduced, tol, max.iter){
	K <- dim(reduced$Mu)[1]
	p <- dim(reduced$Mu)[2]
	n <- dim(X)[1]
	
	la.reduced <- reduced$la
	ind <- as.vector(la.reduced)==0
	
	M <- K - 1 + 2 * K * p + K * p * (p + 1) / 2 - sum(ind)
	reduced.bic <- Manly.bic(reduced$ll, n, M)
	forward.list <- list()
	forward.bic <- rep(NA, sum(ind))
	
	tick <- 0
	for(i in (1:length(ind))[ind]){
		tick <- tick + 1
		la <- la.reduced			for(step in seq(0.01, 0.1, 0.01)){

			la[i] <- step

			res1 <- try(forward <- Manly.EM(X, id = reduced$id, la = la, tol = tol, max.iter = max.iter))
			if(!inherits(res1, "try-error")){

				if(!is.null(forward$ll)){
					break				}

			}
			else {

				la[i] <- -step

				res1 <- try(forward <- Manly.EM(X, id = reduced$id, la = la, tol = tol, max.iter = max.iter))
				if(!inherits(res1, "try-error")){

					if(!is.null(forward$ll)){
						break					}

				}

			}

		}

		if(inherits(res1, "try-error") || (is.null(forward$ll))){
			forward.list[[tick]] <- NULL
			forward.bic[tick] <- Inf
		}
		
		else{	
			forward.list[[tick]] <- forward
			forward.bic[tick] <- Manly.bic(forward$ll, n, M + 1)
		}	

	}

	
	
	return(list(reduced.bic = reduced.bic, forward.bic = forward.bic, forward.list = forward.list))
}



backward <- function(X, full, tol, max.iter){
	K <- dim(full$Mu)[1]
	p <- dim(full$Mu)[2]
	n <- dim(X)[1]
	
	la.full <- full$la
	ind <- as.vector(la.full)!=0
	
	
	M <- K - 1 + K * p + K * p * (p + 1) / 2 + sum(ind)
	full.bic <- Manly.bic(full$ll, n, M)
	reduced.list <- list()
	reduced.bic <- rep(NA, sum(ind))
	
	tick <- 0
	for(i in (1:length(ind))[ind]){
		tick <- tick + 1
		la <- la.full
		la[i] <- 0.0
					
		res1 <- try(reduced <- Manly.EM(X, id = full$id, la = la, tol = tol, max.iter = max.iter))
		if(inherits(res1, "try-error") || (is.null(reduced$ll))){
			reduced.list[[tick]] <- NULL
			reduced.bic[tick] <- Inf
		}
		
		else{	
			reduced.list[[tick]] <- reduced
			reduced.bic[tick] <- Manly.bic(reduced$ll, n, M - 1)
		}
		
	}
	
	return(list(full.bic = full.bic, reduced.bic = reduced.bic, reduced.list = reduced.list))
}

Manly.select <- function(X, model, method, tol = 1e-5, max.iter = 1000, silent = FALSE){

	step <- 0

	if (tol <= 0) stop("Wrong value of tol...\n")
	if (max.iter < 1) stop("Wrong number of iterations iter...\n")


	if(method == "backward"){
		if(!silent){
			repeat{
				step <- step + 1
				C <- backward(X, model, tol, max.iter)
				comparison <- C$reduced.bic < C$full.bic

				if (length(C$reduced.bic != 0)){
					cat("step", step, ":\n\tcurrent BIC =", C$full.bic, "\n\talternative BICs =", C$reduced.bic, "\n")
				} else {
					cat("step", step, ":\n\tcurrent BIC =", C$full.bic, "\n\talternative BICs =", NA, "\n")
				}
				
				if(sum(comparison) != 0){
					index <- which.min(C$reduced.bic)			
					model <- C$reduced.list[[index]]
					BIC <- C$reduced.bic[index]
				} else {
					BIC <- C$full.bic
					break
				}
			}
		}
		else{

			repeat{
				step <- step + 1
				C <- backward(X, model, tol, max.iter)
				comparison <- C$reduced.bic < C$full.bic
				if(sum(comparison)!=0){
					index <- which.min(C$reduced.bic)			
					model <- C$reduced.list[[index]]
					BIC <- C$reduced.bic[index]
				} else {
					BIC <- C$full.bic
					break
				}
			}



		}
	}
	else if(method == "forward"){
		if(!silent){	
			repeat{
				step <- step + 1
				C <- forward(X, model, tol, max.iter)
				comparison <- C$forward.bic < C$reduced.bic 

				if (length(C$forward.bic != 0)){
					cat("step", step, ":\n\tcurrent BIC =", C$reduced.bic, "\n\talternative BICs =", C$forward.bic, "\n")
				} else {
					cat("step", step, ":\n\tcurrent BIC =", C$reduced.bic, "\n\talternative BICs =", NA, "\n")
				}

				if(sum(comparison) != 0){
					index <- which.min(C$forward.bic)			
					model <- C$forward.list[[index]]
					BIC <- C$forward.bic[index]
				} else {
					BIC <- C$reduced.bic
					break
				}


			}
		}
		else{

			repeat{
				step <- step + 1
				C <- forward(X, model, tol, max.iter)
				comparison <- C$forward.bic < C$reduced.bic 
				if(sum(comparison) != 0){
					index <- which.min(C$forward.bic)			
					model <- C$forward.list[[index]]
					BIC <- C$forward.bic[index]
				} else {
					BIC <- C$reduced.bic
					break
				}


			}
		}
		
	}
	return(list(la = model$la, tau = model$tau, Mu = model$Mu, S = model$S, gamma = model$gamma, id = model$id, ll = model$ll, bic = BIC, iter = model$iter, flag = model$flag))
}



