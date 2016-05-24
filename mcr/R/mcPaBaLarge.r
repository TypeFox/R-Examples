# TODO: Interface to C-implementation of the Passing-Bablok algorithm for large datasets.
#       This adaption of the Passing-Bablok algorithm was implemented in DLLs formerly
#       used for method comparison procedures in winCAEv, written by Dr. Carl Kuhn (Roche
#       Diagnostics, Penzberg, Department Biostatistics & Online Data Processing).
#       This C++ source-code is now reimplemented in C and adapted to gain more flexibility. 
#       
# Author: schueta6
###############################################################################



#' Passing-Bablok Regression for Large Datasets
#' 
#' This function represents an interface to a fast C-implementation of an adaption of the Passing-Bablok
#' algorithm for large datasets. Instead of building the complete matrix of pair-wise slope values, a pre-defined
#' binning of slope-values is used (Default NBins=1e06). This reduces the required memory dramatically and speeds
#' up the computation.
#'
#' @param X         (numeric) vector containing measurement values of reference method
#' @param Y         (numeric) vector containing measurement values of test method
#' @param NBins     (integer) value specifying the number of bins used to classify slope-values
#' @param alpha     (numeric) value specifying the 100(1-alpha)\% confidence level for confidence intervals
#' @param posCor    (logical) should algorithm assume positive correlation, i.e. symmetry around slope 1?
#' @param calcCI    (logical) should confidence intervals be computed?
#' 
#' @return Matrix of estimates and confidence intervals for intercept and slope. No standard errors provided by this algorithm. 
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com} (partly re-using code of function 'mc.paba')
#' 
#' @examples
#' 
#'  library("mcr")
#'  data(creatinine,package="mcr")
#'  
#' # remove any NAs
#' crea <- na.omit(creatinine)
#'     
#' # call the approximative Passing-Bablok algorithm (Default NBins=1e06)
#' res1 <- mcreg(x=crea[,1], y=crea[,2], method.reg="PaBaLarge", method.ci="analytical") 
#' getCoefficients(res1)
#' 
#' # now increase the number of bins and see whether this makes a difference
#' res2 <- mcreg(x=crea[,1], y=crea[,2], method.reg="PaBaLarge", method.ci="analytical", NBins=1e07) 
#' getCoefficients(res2)
#' getCoefficients(res1)-getCoefficients(res2)


mc.paba.LargeData <- function(X, Y, NBins=1e06, alpha=0.05, posCor=TRUE, calcCI=TRUE) 
{
    NBinsUpperBound <- 1e08                     # NBins may not be larger than that
    
    ## Check validity of parameters
    
    stopifnot(length(X)==length(Y))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(!is.na(X))
    stopifnot(!is.na(Y))
    stopifnot(!is.na(posCor))
    
    alpha <- as.numeric(alpha)
    stopifnot(!is.na(alpha))
    stopifnot( alpha > 0 && alpha < 1)
    
    NBins <- as.integer(NBins)
    stopifnot(!is.na(NBins))
        
    if(NBins > NBinsUpperBound)
        stop("The maximum number of bins is '", NBinsUpperBound,"' and may not be exceeded!")
    
    ## define double constants which serves as surrogate value for Inf and -Inf which are not allowed in calls to C-functions
    
    INF  <-  1020304.050607                    # this value should be unique enough, such that no computed values will ever be equal to it
    
    # Call the C-function 
    
    res <- .C( "PaBaLargeData", 
                pX=as.double(X), pY=as.double(Y), pNData=as.integer(length(X)),                                     # input
                pPosCor=as.integer(posCor), pNBins=as.integer(NBins),
                pSlope=double(1), pQuantile=as.double(qnorm(1-alpha/2)),                                            # output
                pSlopeLower=as.double(-INF), pSlopeUpper=as.double(INF), pCIundefined=as.integer(0))

    mcres.slope <- res$pSlope  

## Confidence intervals for slope (only keep results coming from the C-function if they are not equal to -INF and INF)

    if(calcCI) 
    {
        if(res$pCIundefined == 1)                   # errors occurred -> CI for slope is undefined -> CI for intercept also undefined
        {
            mcres.slopeL <- as.numeric(NA)
            mcres.slopeU <- as.numeric(NA)
        }
        else
        {
            mcres.slopeL <- ifelse(res$pSlopeLower == -INF, -Inf, res$pSlopeLower)
            mcres.slopeU <- ifelse(res$pSlopeUpper == INF, Inf, res$pSlopeUpper)
        }
    }
    else 
    {
        ## No theoretical CIs computed, use resampling to get CIs
        mcres.slopeL <- as.numeric(NA)
        mcres.slopeU <- as.numeric(NA)
    }    

    ## Intercept

    mcres.intercept <- median(calcDiff(Y, mcres.slope*X))
    if(calcCI) 
    {
        if(res$pCIundefined == 1)                   # errors occurred -> CI for slope is undefined -> CI for intercept also undefined
        {
            mcres.interceptL <- as.numeric(NA)
            mcres.interceptU <- as.numeric(NA)
        }
        else
        {
            if(mcres.slopeL==-Inf) 
                mcres.interceptL <- -Inf
            else 
                mcres.interceptL <- median(calcDiff(Y,mcres.slopeU*X))
            
            if(mcres.slopeU==Inf) 
                mcres.interceptU <- Inf
            else 
                mcres.interceptU <- median(calcDiff(Y,mcres.slopeL*X))
        }
    }
    else 
    {
        mcres.interceptL <- as.numeric(NA)
        mcres.interceptU <- as.numeric(NA)
    }
    
    ## Prepare result matrix
    rmat <- matrix(nrow=2,ncol=4)
    rownames(rmat) <- c("Intercept","Slope")
    colnames(rmat) <- c("EST","SE","LCI","UCI")
    rmat[,1] <- c(mcres.intercept,mcres.slope)
    rmat[,2] <- NA
    rmat["Intercept","LCI"] <- mcres.interceptL
    rmat["Intercept","UCI"] <- mcres.interceptU
    rmat["Slope","LCI"] <- mcres.slopeL
    rmat["Slope","UCI"] <- mcres.slopeU

    return(rmat)
}