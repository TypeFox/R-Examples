# Provides two methods to binarize a vector consisting of real values. The methods are scale-space-based. With <method>="A", 
# the threshold for the binarization is calculated according to the BASC A algorithm and with <method>="B" according to the 
# BASC B algorithm. For details on how this is done have a look at the vignette. The <tau> argument is a parameter for the bootstrap-test
# and it is an indicator what quality the binarization should have. The resulting p-value is shows how good the binarization fulfills this
# quality requirement. With <numberofSamples> you can control the number of samples that are used for the bootstrap test.
# <sigma> is only used for the BASC B algorithm and so it will be ignored if <method> equals "A". If <method>="B" then <sigma> should be a sequence
# of ascending numbers and the values are used as parameters for the bessel function.
binarize.BASC <- function(vect, method=c("A","B"), tau = 0.01, numberOfSamples = 999, sigma = seq(0.1, 20, by=.1), na.rm=FALSE){    
    # Check type of input vector
    if(!is.numeric(vect))
        stop("The input vector must consist of numerical values!")    
    
    # Check for NA values
    if (!na.rm && any(is.na(vect)))
      stop("Cannot binarize in the presence of NA values!")
    else
    if (na.rm)
    {
      vect <- vect[!is.na(vect)]
    }
    
    if (any(!is.finite(vect)))
     stop("Cannot binarize Inf values!")

    
    # Check type of method
    method <- match.arg(method, c("A","B"))
    
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
        res <- binarizeBASCA_C(vect, tau, numberOfSamples)
        
        result <- new(
                    "BASCResult",
                    originalMeasurements = vect, 
                    binarizedMeasurements = res$binarized_vector, 
                    threshold = res$threshold, 
                    p.value = res$p_value,
                    intermediateSteps = res$other_results$P_Mat,
                    intermediateHeights = res$other_results$H_Mat,
                    intermediateStrongestSteps = res$other_results$v_vec,
                    method = "BASC A"
                )
    }
    else if(method == "B"){
        res <- binarizeBASCB_C(vect, tau, numberOfSamples, sigma)
        
        result <- new(
                        "BASCResult",
                        originalMeasurements = vect,  
                        binarizedMeasurements = res$binarized_vector, 
                        threshold = res$threshold, 
                        p.value = res$p_value,
                        intermediateSteps = res$other_results$steps,
                        intermediateHeights = res$other_results$H_Mat,
                        intermediateStrongestSteps = res$other_results$v_vec,
                        method = "BASC B"
                    )
    }
    else{
        stop(sprintf("'method' has to be either \"A\" or \"B\"!", method))
    }    
    
    return(result)
}

#interface for the Call to the C-Function, which does all the calculations for binarization.BASC with method
#argument = "A" 
binarizeBASCA_C <- function(vect, tau, numberOfSamples){        
    
    #call the C-Function
    result <- .Call("binarizeBASCA", 
                    as.double(vect), 
                    as.double(tau), 
                    as.integer(numberOfSamples))
    
    #all the matrices have to be transposed, due to different interpretation of the sequence of values
    #between C and R (in C line by line, and in R column by column)
    result$other_results$Cc <- t(result$other_results$Cc)
    result$other_results$Ind <- t(result$other_results$Ind)
    result$other_results$P_Mat <- t(result$other_results$P_Mat)
    result$other_results$Q_Mat <- t(result$other_results$Q_Mat)
    result$other_results$H_Mat <- t(result$other_results$H_Mat)
    
    return(result);
}

#interface for the Call to the C-Function, which does all the calculations for binarization.BASC with method
#argument = "B"
binarizeBASCB_C <- function(vect, tau, numberOfSamples, sigma){
    
    result <- .Call("binarizeBASCB", 
                    as.double(vect), 
                    as.double(tau), 
                    as.integer(numberOfSamples), 
                    as.double(sigma))
    
    #all the matrices have to be transposed, due to different interpretation of the sequence of values
    #between C and R (in C line by line, and in R column by column)
    result$other_results$smoothed = t(result$other_results$smoothed)
    result$other_results$zerocrossing = t(result$other_results$zerocrossing)
    result$other_results$steps = t(result$other_results$steps)
    result$other_results$H_Mat = t(result$other_results$H_Mat)
    result$other_results$smoothedX = t(result$other_results$smoothedX)
    result$other_results$meanlist = t(result$other_results$meanlist)
    
    return(result)
}
