TASC <- function(vect,
		method=c("A","B"),
		tau = 0.01,
		numberOfSamples = 999,
		sigma = seq(0.1, 20, by=.1),
		na.rm=FALSE,
		error = c("mean", "min")){    
	
	# Check type of input vector
	if(!is.numeric(vect))
		stop("The input vector must consist of numerical values!")    
	
	# Check for NA values
	if (!na.rm && any(is.na(vect)))
	  stop("Cannot trinarize in the presence of NA values!")
	else
	if (na.rm)
	{
	  vect <- vect[!is.na(vect)]
	}
	
	if (any(!is.finite(vect)))
	 stop("Cannot trinarize Inf values!")

	# Check type of method
	method <- match.arg(method, c("A","B"))
	error <- match.arg(error, c("mean", "min"))
	
	# Check type and value of tau
	if(!is.numeric(tau))
		stop("'tau' must be numeric!")
	if(tau < 0 || tau > 1)
		stop("'tau' has to be in [0,1]!")
		
	# Check type and value of numberofSamples
	if(!is.numeric(numberOfSamples))
		stop("'numberOfSamples' must be numeric!")
	if(numberOfSamples < 0)
		stop("'numberOfSamples' has to be >= 0!")
		
	# Check type of sigma
	if(method == "B" && !is.numeric(sigma))
		stop("'sigma' must consist of numerical values!")
	
	# Verify input length
	if(length(vect) < 3)
		stop("The input vector must have at least 3 entries!")

	# Check whether the input is constant
	if(length(unique(vect)) == 1)
		stop("The input vector is constant!")
		
		
		
	#set.seed(proc.time()[3]*1000)
	runif(1)
	
 
	#calculate the results according to the <method> argument
	if(method == "A"){
		if(error == "min"){
			res <- TASCA_C_min(vect, tau, numberOfSamples)
			me <- "TASC A (min)"
		}else{
			res <- TASCA_C(vect, tau, numberOfSamples)
			me <- "TASC A"
		}
		
		result <- new(
					"TASCResult",
					originalMeasurements = vect, 
					trinarizedMeasurements = res$binarized_vector, 
					threshold1 = res$threshold1,
					threshold2 = res$threshold2, 
					p.value = res$p_value,
					intermediateSteps = res$other_results$P_Mat[-1,,drop=FALSE],
					intermediateHeights1 = res$other_results$H_Mat1[-1,,drop=FALSE],
					intermediateHeights2 = res$other_results$H_Mat2[-1,,drop=FALSE],
					intermediateStrongestSteps = res$other_results$v_vec[-1,,drop=FALSE],
					method = me
				)
	}else if(method == "B"){
		if(error == "min"){
			res <- TASCB_C_min(vect, tau, numberOfSamples, sigma)
			me <- "TASC B (min)"
		}else{
			res <- TASCB_C(vect, tau, numberOfSamples, sigma)
			me <- "TASC B"
		}
		
		result <- new(
					"TASCResult",
					originalMeasurements = vect, 
					trinarizedMeasurements = res$binarized_vector, 
					threshold1 = res$threshold1,
					threshold2 = res$threshold2, 
					p.value = res$p_value,
					intermediateSteps = res$other_results$steps[-1,,drop=FALSE],
					intermediateHeights1 = res$other_results$H_Mat1[-1,,drop=FALSE],
					intermediateHeights2 = res$other_results$H_Mat2[-1,,drop=FALSE],
					intermediateStrongestSteps = res$other_results$v_vec[-1,,drop=FALSE],
					method = me
				)
	}else{
		stop(sprintf("'method' has to be either \"A\" or \"B\"!", method))
	}
	
	return(result)
}

TASCA_C <- function(vect, tau, numberOfSamples){
	#call the C-Function
	result <- .Call("TASCA", 
					as.double(vect), 
					as.double(tau), 
					as.integer(numberOfSamples))
	
	#all the matrices have to be transposed, due to different interpretation of the sequence of values
	#between C and R (in C line by line, and in R column by column)
	result$other_results$Cc <- t(result$other_results$Cc)
	result$other_results$Ind <- t(result$other_results$Ind)
	result$other_results$P_Mat <- t(result$other_results$P_Mat)
	result$other_results$Q_Mat <- t(result$other_results$Q_Mat)
	result$other_results$H_Mat1 <- t(result$other_results$H_Mat1)
	result$other_results$H_Mat2 <- t(result$other_results$H_Mat2)
	result$other_results$v_vec <- matrix(result$other_results$v_vec, ncol = 2, byrow = TRUE)

	return(result);
}


TASCB_C <- function(vect, tau, numberOfSamples, sigma){
    
    result <- .Call("TASCB", 
                    as.double(vect), 
                    as.double(tau), 
                    as.integer(numberOfSamples), 
                    as.double(sigma))
    
    #all the matrices have to be transposed, due to different interpretation of the sequence of values
    #between C and R (in C line by line, and in R column by column)
    result$other_results$smoothed = t(result$other_results$smoothed)
    result$other_results$zerocrossing = t(result$other_results$zerocrossing)
    result$other_results$steps = t(result$other_results$steps)
    result$other_results$H_Mat1 = t(result$other_results$H_Mat1)
    result$other_results$H_Mat2 = t(result$other_results$H_Mat2)
    result$other_results$smoothedX = t(result$other_results$smoothedX)
    result$other_results$meanlist = t(result$other_results$meanlist)
    result$other_results$v_vec <- matrix(result$other_results$v_vec, ncol = 2, byrow = TRUE)
    
    return(result)
}

TASCA_C_min <- function(vect, tau, numberOfSamples){
	#call the C-Function
	result <- .Call("TASCA_min", 
					as.double(vect), 
					as.double(tau), 
					as.integer(numberOfSamples))
	
	#all the matrices have to be transposed, due to different interpretation of the sequence of values
	#between C and R (in C line by line, and in R column by column)
	result$other_results$Cc <- t(result$other_results$Cc)
	result$other_results$Ind <- t(result$other_results$Ind)
	result$other_results$P_Mat <- t(result$other_results$P_Mat)
	result$other_results$Q_Mat <- t(result$other_results$Q_Mat)
	result$other_results$H_Mat1 <- t(result$other_results$H_Mat1)
	result$other_results$H_Mat2 <- t(result$other_results$H_Mat2)
	result$other_results$v_vec <- matrix(result$other_results$v_vec, ncol = 2, byrow = TRUE)

	return(result);
}


TASCB_C_min <- function(vect, tau, numberOfSamples, sigma){
    
    result <- .Call("TASCB_min", 
                    as.double(vect), 
                    as.double(tau), 
                    as.integer(numberOfSamples), 
                    as.double(sigma))
    
    #all the matrices have to be transposed, due to different interpretation of the sequence of values
    #between C and R (in C line by line, and in R column by column)
    result$other_results$smoothed = t(result$other_results$smoothed)
    result$other_results$zerocrossing = t(result$other_results$zerocrossing)
    result$other_results$steps = t(result$other_results$steps)
    result$other_results$H_Mat1 = t(result$other_results$H_Mat1)
    result$other_results$H_Mat2 = t(result$other_results$H_Mat2)
    result$other_results$smoothedX = t(result$other_results$smoothedX)
    result$other_results$meanlist = t(result$other_results$meanlist)
    result$other_results$v_vec <- matrix(result$other_results$v_vec, ncol = 2, byrow = TRUE)
    
    return(result)
}
