"lnormal.ssf" <- function(method = c("1", "2", "3", "4"), fixed, loge, useFixed = FALSE)
{
    method <- match.arg(method)
    
    ## Defining helper functions (used below)
    ytrans <- function(y, cVal, dVal) {qnorm((y-cVal)/(dVal-cVal))}
    bfct <- function(x, y, cVal, dVal, eVal) {ytrans(y, cVal, dVal)/log(x/eVal)}
    efct <- function(x, y, bVal, cVal, dVal) {x * exp(-ytrans(y, cVal, dVal)/bVal)}
    
    ## Assigning function for finding initial b and e parameter values    
    findbe <- switch(method,
    "1" = findbe1(function(x) {rVec <- log(x); rVec[!x>0] <- NA; rVec}, 
    function(y, cVal, dVal) {qnorm((dVal - y) / (1.01 * max(dVal - y)))}, -1),    
    "2" = findbe2(bfct, efct, "Anke", -1),
    "3" = findbe3(-1),
    "4" = findbe2(bfct, efct, "Normolle", -1))
    
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        ## Finding initial values for c and d parameters
        cdVal <- findcd(x, y)
        if (useFixed) {}  # not implemented at the moment        
    
        ## Finding initial values for b and e parameters    
        beVal <- findbe(x, y, cdVal[1], cdVal[2])       
        if (loge) {beVal[2] <- log(beVal[2])}

        return(c(beVal[1], cdVal, beVal[2])[is.na(fixed)])
    }
}