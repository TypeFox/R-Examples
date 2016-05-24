################################################################################
# TODO LIST
# TODO: ...

################################################################################
# NOTES
# Stutter simulation require package  'mc2d' + dep. 'mvtnorm'

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 28.10.2014: First version.

#' @title Vectorized Multinomial Distribution For Large Numbers
#'
#' @description
#' Simulates one cycle of the PCR process.
#'
#' @details
#' Generate multinomially distributed random number vectors using \code{rmultinomial}.
#' Handles integer overflow by splitting numbers in chunks and call \code{rmultinomial} repeatedly.
#' 
#' @param n numeric vector with number of observations.
#' @param size numeric vector specifying the number of molecules going into the PCR cycle.
#' @param prob numeric non-negative matrix of size (x x K) specifying the probability for the K classes.
#' @param check logical for extended checking of parameters.
#' @param debug logical to print debug information.
#' 
#' @return matrix with simulated results
#' 
#' @importFrom mc2d rmultinomial
#' @importFrom utils head
#' 

rmultinomxl <- function(n, size, prob, check=TRUE, debug=FALSE){
  
  # Debug info.
  if(debug){
    print(paste(">>>>>> IN:", match.call()[[1]]))
    print("CALL:")
    print(match.call())
    print("###### PROVIDED ARGUMENTS")
    print("n:")
    print(head(n))
    print("size:")
    print(head(size))
    print("prob:")
    print(head(prob))
  }
  
  # CHECK PARAMETERS ##########################################################

  # Extended check.
  if(check){

    if(n != length(size)){
      stop(paste("'n' must be equal to length of 'size'."))
    }
    
  }
  
  # Get max integer size.
  .imax <- .Machine$integer.max
  
  # To avoid error caused by integer overflow:
  # Divide in max integer chunks and add up after rmultinomial.
  iChunks <- floor(size / .imax)
  
  if(debug){
    print("iChunks")
    print(head(iChunks))
  }
  
  if(all(iChunks == 0)) {
    # If number of molecules is < .imax.

    if(debug){
      print(paste("size <", .imax, "Running rmultinomial one time."), sep="")
    }
    
    res <- mc2d::rmultinomial(n=n, size=size, prob=prob)
    
  } else {
    # If number of .imax chunks is at least one...

    # ... run rmultinomial for all chunks.
    maxchunks <- max(iChunks)
    
    if(debug){
      print(paste("size >", .imax, ". Forced to run rmultinomial ", maxchunks, " times."), sep="")
    }

    # Calculate the rest. 
    rest <- size - iChunks * .imax
    rest[rest < 0] <- 0
    
    if(debug){
      print("rest")
      print(head(rest))
    }
    
    # Then add the rest (for #molecules < .imax this is the only run).
    res <- mc2d::rmultinomial(n=n, size=rest, prob=prob)

    if(debug){
      print(paste("Result for the rest:"))
      print(head(res))
    }
    
    for(i in 1:maxchunks){

      if(debug){
        print(paste("Chunk", i, "of", maxchunks))
      }
      
      # Mask places which does not overflow (and hence should not get more molecules).
      mask <- iChunks
      mask[iChunks <= 0] <- 0
      mask[iChunks > 0] <- 1

      # Draw molecules with the given probabilities.
      tmp <- mc2d::rmultinomial(n=n, size=.imax, prob=prob)
      
      if(debug){
        print("tmp")
        print(head(tmp))
      }
      
      # Add to result.
      res <- res + (tmp * mask)
 
      if(debug){
        print(paste("Mask vector:", head(paste(mask, collapse=","))))
        print(paste("Result:"))
        print(head(res))
      }
      
      # Reduce iChunks (to update mask).
      iChunks <- iChunks - 1
      
    }
   
  }
  
  # RETURN ####################################################################
  
  # Debug info.
  if(debug){
    print(paste("<<<<<< EXIT:", match.call()[[1]]))
  }
  
  return(res)
  
}