# global enviroment:
rrecsys.env <- new.env()
rrecsys.env$nrLoops <- 10
rrecsys.env$autoConverge <- FALSE
rrecsys.env$counter <- 0
rrecsys.env$last_error <- 0
rrecsys.env$deltaErrorThreshold <- 1e-05
rrecsys.env$showError <- FALSE
rrecsys.env$minNrLoops <- 10

setStoppingCriteria <- function(autoConverge = FALSE, deltaErrorThreshold = 1e-05, nrLoops = NULL, minNrLoops = 10) {
    if (!autoConverge & is.null(nrLoops)) {
        cat("Please specify either autoConvergence or nrLoops criteria.")
        return(showStoppingCriteria())
    }
    if (autoConverge) {
        
        rrecsys.env$autoConverge <- autoConverge
        rrecsys.env$deltaErrorThreshold <- deltaErrorThreshold
        rrecsys.env$minNrLoops <- minNrLoops
        if (minNrLoops < 1) 
            stop("Invalid number of minNrLoops!")
        return(showStoppingCriteria())
    }
    if (missing(autoConverge) & !missing(nrLoops)) {
        if (nrLoops < 1) 
            stop("Invalid number of loops!")
        rrecsys.env$nrLoops <- nrLoops
        rrecsys.env$autoConverge <- FALSE
        return(showStoppingCriteria())
    }
}

showStoppingCriteria <- function() {
    
    if (!rrecsys.env$autoConverge) {
        cat("Stopping criteria is defined by the number of loops. Current loop number is: ", rrecsys.env$nrLoops, ".")
    } else {
        cat("Stopping criteria is defined by a threshold on the global delta error(in terms of difference of RMSE on two consecutive iterations).\nWARNING: Current configuration can diverge.You can show or not delta error(in terms of difference of global RMSE on two consecutive iterations) over time by toggling on/off method showError().\nThe threshold is: ", 
            rrecsys.env$deltaErrorThreshold, ".")
    }
    
    cat("\nAlgorithm that support autoConvergence: FunkSVD, wALS, BPR.")
    
    if (rrecsys.env$showError) {
        cat("\nShow delta error is ON!")
    } else {
        cat("\nShow delta error is OFF!")
    }
    
}

showDeltaError <- function() {
    
    if (rrecsys.env$showError) {
        rrecsys.env$showError <- FALSE
        cat("Show delta error is OFF!")
    } else {
        rrecsys.env$showError <- TRUE
        cat("Show delta error is ON!")
    }
}


resetrrecsysenv <- function() {
    rrecsys.env$counter <- 0
    rrecsys.env$last_error <- 0
}


isConverged <- function(x, p) {
    
  #in terms of RMSE
  error <- sqrt(sum(abs(x - p)^2)/length(x))
  delta_error <- abs(rrecsys.env$last_error - error)

  if(is.infinite(delta_error) || is.nan(delta_error)) stop("Error diverges!!! Fix learning rate or regularization term.")
  
    if (rrecsys.env$autoConverge) {
        # error calculated in terms of RMSE
        
        rrecsys.env$counter <- rrecsys.env$counter + 1
        if ((rrecsys.env$counter > rrecsys.env$minNrLoops) & (delta_error < rrecsys.env$deltaErrorThreshold)) {
            return(TRUE)
        }

        if (rrecsys.env$showError) 
            writeLines(sprintf("Iteration: %s. Delta error: %s.", rrecsys.env$counter, delta_error))
        
        # supose there are 10 consecutive iteration where the delta error never changes, we consider our iteration converged if (rrecsys.env$counter == 10) areWeDone <-
        # TRUE
        
        rrecsys.env$last_error <- error
    } else {
        
        if (rrecsys.env$showError) {
            
            rrecsys.env$last_error <- error
            writeLines(sprintf("Iteration: %s. Delta error: %s.", rrecsys.env$counter, delta_error))
        }
        
        
        rrecsys.env$counter <- rrecsys.env$counter + 1
        if (rrecsys.env$counter == rrecsys.env$nrLoops) 
            return(TRUE)
        
    }
  
    return(FALSE)
} 
