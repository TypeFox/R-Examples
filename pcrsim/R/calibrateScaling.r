################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 30.10.2015: First version.

#' @title Calibrate Peak Height Scaling
#'
#' @description
#' Corrects the peak height scaling intercept.
#'
#' @details
#' Finds the scaling factor correction that gives the best fit between the
#' simulated peak heights and the target average peak height.
#' This is the final step in calibrating pcrsim for a specific kit and method.
#' 
#' @param data data.frame simulated data.
#' @param ref data.frame with reference profiles for the simulated data.
#'  If NULL the best guess will be used (see \code{\link{guessProfile}}).
#' @param target numeric target average peak height 'H'.
#' @param c.min numeric the smallest correction factor in the current range.
#' @param c.max numeric the largest correction factor in the current range.
#' @param min.step.size numeric threshold. Exit function when step size is this small.
#' @param exact.matching logical to indicate exact reference sample to dataset sample name matching.
#' @param progress logical flag to show progress messages.
#' @param .i integer internal counter for number of iterations.
#' 
#' @return numeric corrected scaling intercept.
#' 
#' @importFrom strvalidator filterProfile calculateHeight addData
#' 
#' @export
#' 

calibrateScaling <- function(data, ref=NULL, target, c.min=0, c.max=10000,
                             min.step.size=0.0001, exact.matching = FALSE,
                             progress=FALSE, .i=1) {
  
  # CHECK PARAMETERS ##########################################################
  
  if(!is.data.frame(data)){
    stop(paste("'data' must be a data.frame"))
  }
  
  if(!all(c("Height", "CE.RFU") %in% names(data))){
    stop(paste("'data' must contain a column 'Height' and a column 'CE.RFU'"))
  }
  
  if(!any(is.data.frame(ref), is.null(ref))){
    stop(paste("'ref' must be a data.frame or NULL"))
  }
  
  if(!all(is.numeric(target), target > 0)){
    stop(paste("'target' must be a positive numeric"))
  }
  
  if(!all(is.numeric(c.min), c.min >= 0)){
    stop(paste("'c.min' must be a positive numeric >= 0"))
  }
  
  if(!all(is.numeric(c.max), c.max > c.min)){
    stop(paste("'c.max' must be a positive numeric > 'c.min'"))
  }
  
  if(!all(is.numeric(min.step.size), min.step.size > 0)){
    stop(paste("'min.step.size' must be a positive numeric"))
  }

  if(!is.logical(progress)){
    stop(paste("'progress' must be logical."))
  }
  
  # PREPARE ###################################################################

  # Declare variables.
  lastrss <- Inf  # Last rss value, initated to Inf so that all new rss < lastrss.
  rssfit <- NA    # rss for best fit (so far).
  maxit <- 10      # Max iterations.
  finished <- FALSE # Flag end of run.

  decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(x+1)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  
  # Get number of decimal places.
  dec <- decimalplaces(min.step.size)

  # If first iteration perform some tasks.
  if(.i == 1){

    # Check if reference is provided.
    if(is.null(ref)){
      
      # Guess the correct profile.
      ref <- strvalidator::guessProfile(data = data, height = 50, ratio = 0.5,
                                        na.rm = TRUE, ol.rm = TRUE)
      
      # Use exact matching when adding information.
      exact.matching <- TRUE
      message("exact.matching set to TRUE since ref=NULL")
      
    }

    # Filter data.
    data <- strvalidator::filterProfile(data=data, ref=ref, exact=exact.matching)
    
    # Add heterozygous indicator.
    ref <- strvalidator::calculateHeterozygous(ref)
    data <- strvalidator::addData(data=data, new.data = ref,
                                  by.col = "Sample.Name", then.by.col = "Marker",
                                  exact = exact.matching)

    if(progress){
      print("Please review the data to make sure all profiles have been correctly filtered:")  
      print(strvalidator::calculateHeight(data = data, add = FALSE))
    }
    
  }
  
  # OPTIMISE ##################################################################
  
  # Print progress.
  if(progress){
    message(paste("\nIteration ", .i,
                  ", optimise in range ", c.min, "-", c.max, ":", sep=""))
  }
  
  # Calculate correction factor step size.  
  cStep <- (c.max - c.min) / 50
  
  # Create a uniformly distributed sequence.
  cVector <- seq(c.min, c.max, cStep)
  cFit <- NA # Initiate to NA.
  
  # Repeat over the sequence.
  for(d in seq(along=cVector)){
    
    # Adjust height using the factor.
    data$Height <- data$CE.RFU * cVector[d]

    # Calculate average peak height.
    dataH <- strvalidator::calculateHeight(data = data, add = FALSE)
    
    # Calculate the residual sum of squares.
    rss <- sum((dataH$H - target)^2)
    
    # Check fit.
    if(!is.na(rss) & !is.na(lastrss)){
      
      if (rss < lastrss){
        
        # Save current values.
        cFit <- cVector[d]
        rssfit <- rss
        
        # Calculate overall corrected mean average peak height.
        meanH <- mean(dataH$H)
        
        # Print progress.
        if(progress){
          message(paste("Current correction factor = ", round(cFit, dec),
                        ", RSS = ", round(rss),
                        ", mean(H) = ", round(meanH), sep=""))
        }
        
      }
    }
    
    # Save current rss value.
    lastrss <- rss
    
  } # End of for loop.
  
  # Recursive if not optimised enough.
  if(.i == maxit){
    
    # Print progress.
    if(progress){
      message(paste("\nMaximum number of iterations (", maxit ,") reached!", sep=""))
      message(paste("Best fit correction factor = ", cFit, sep=""))
    }
    
  } else if(!is.na(rssfit) && cStep > min.step.size){
    
    # Calculate new range and fix to {0,}.
    cNewMin <- cFit - cStep
    cNewMax <- cFit + cStep
    if(cNewMin < 0){
      cNewMin <- 0
    }

    # Optimise further.
    cFit <- calibrateScaling(data=data, ref=NULL, target=target,
                             c.min=cNewMin, c.max=cNewMax,
                             min.step.size=min.step.size, progress=progress,
                             .i=.i+1)
    
  } else if(!is.na(rssfit) && cStep < min.step.size){
    
    # Print progress.
    if(progress){
      message(paste("\nMinimum step size (", min.step.size, ") reached!", sep=""))
      message(paste("Best fit correction factor = ", cFit, sep=""))
    }
    
    
  } else {
    
    # This should newer happen?
    message(paste("This was unexpected! Dumping some data:", sep=""))
    message(paste("optimise in range ", c.min, "-", c.max, ":", sep=""))
    message(paste("iteration ", .i, sep=""))
    message(paste("rssfit ", rssfit, sep=""))
    message(paste("Reached step size = ", cStep, sep=""))
    message(paste("Best fit correction factor = ", cFit, sep=""))

  }
  
  if(.i==1){
    # This is the last thing happening in the recursive function.

    # Get original intercept.    
    intercept <- unique(data$CE.S.Intercept)
    
    # Calculate corrected scaling intercept. 
    cFit <- intercept + log(cFit)

    # Print progress.
    if(progress){
      message(paste("Corrected scaling intercept = ", cFit, sep=""))
    }

  }

  return (cFit)
  
}

