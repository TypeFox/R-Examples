# File wquantile.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2010-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                    'wquantile'                                               #
################################################################################                                                                                
# Purpose: This function computes weighted quantiles of each column (by default# 
#          or for each row if specified by the user) of a matrix/data.frame    #
#          It is a wrapper to the 'wtd.quantiles' function of the 'Hmisc'      #
#          package, specially thought for a matrix containing streamflows      #
#          simulated by  different (behavioural) parameter sets                #                             
################################################################################                                                                                
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: June 09th, 2010                                                     #        
# Updates: 27-Sep-2013                                                         #
################################################################################ 
# 

# Result\Value: It returns a matrix with an amount of rows equal to the
# amount of columns in 'x', and as many columns (by default, or for each 
# row if specified by the user) as 'length(probs)'

# Original idea taken from: http://rwiki.sciviews.org/doku.php?id=guides:tutorials:hydrological_data_analysis:glue

# x       : a numeric matrix, which columns (by default) contains the values to be 
#           used for the computation of the weighted quantiles
# weights : a numeric vector of weights. \cr
#           Omitting the 'weights' argument or specifying 'NULL' or a 
#           zero-length vector will result in the usual unweighted estimates.
# byrow   : logical, indicating if the computations have to be made for each column or for each row of \code{x}.
#           When the values to be used for the computation of the quantiles are stored in rows, \code{byrow} must be \kbd{TRUE}. \cr
#           When the values to be used for the computation of the quantiles are stored in columns, \code{byrow} must be \kbd{FALSE}.
# probs   : a vector of quantiles to compute.  
#           Default is c(.025, .5, .975) (2.5\%, 50\%, 97.5\%)
# normwt  : specify 'normwt=TRUE' to make 'weights' sum to 'length(x)' after deletion of NAs
# verbose : logical; if TRUE, progress messages are printed

wquantile <- function(x, weights=NULL, byrow=FALSE, probs=c(.025, .5, .975), 
                      normwt=TRUE, verbose=TRUE ) {
  
  # computation of the 95PPU  of the Parameters
  #require(Hmisc) #   'wtd.quantile'
  
  if (byrow==TRUE) {
  
    nobs <- ncol(x)
    n    <- nrow(x)
    
  } else { 
      
      nobs <- nrow(x) 
      n    <- ncol(x)
    
    } # ELSE end
  
  # If the user provided the weights
  if ( !missing(weights) ) {
   if ( length(weights) != nobs ) {
       if (byrow==TRUE) {
         stop( paste("Invalid argument: 'length(weights) != ncol(x)' (", length(weights), "!=", nobs, ")", sep="" ) ) 
       } else stop( paste("Invalid argument: 'length(weights) != nrow(x)' (", length(weights), "!=", nobs, ")", sep="" ) ) 
   } # IF end
  } # IF end
  
  if (verbose) pbar <- txtProgressBar(min = 0, max = n, style = 3)
  
  tmp <- sapply(1:n, function(i, p) {
  
         #if (verbose) print( paste(i, "/", n, format(round(100*i/n, 2), width=6, justify="left"), 
         #                             "%", sep=""), quote=FALSE )     
         if (verbose) setTxtProgressBar(pbar, i)
         
         if (byrow==TRUE) {
         
           as.numeric(Hmisc::wtd.quantile(p[i, ], weights = weights, probs = probs, normwt=normwt))
           
         } else  as.numeric(Hmisc::wtd.quantile(p[, i], weights = weights, probs = probs, normwt=normwt))

         }, p = x)
         
  # clossing the progress bar
  if (verbose) close(pbar)
  
  p95 <- t(tmp)
  colnames(p95) <- paste(probs*100, "%", sep="")  
  
  if (byrow==TRUE) {
    if ( !is.null(rownames(x)) ) rownames(p95) <- rownames(x)
  } else {
        if ( !is.null(colnames(x)) )  rownames(p95) <- colnames(x)
    } # ELSE end
  
  return(p95)  

} #'wquantile' END
