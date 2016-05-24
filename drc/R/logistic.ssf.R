"logistic.ssf" <- function(method = c("1", "2", "3", "4"), fixed, useFixed = FALSE)
{
    method <- match.arg(method)
    
    ## Defining helper functions (used below)
    ytrans <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
    bfct <- function(x, y, cVal, dVal, eVal) {ytrans(y, cVal, dVal) / (x - eVal)}
    efct <- function(x, y, bVal, cVal, dVal) {x - ytrans(y, cVal, dVal) / bVal}
    
    ## Assigning function for finding initial b and e parameter values    
    findbe <- switch(method,
    "1" = findbe1(function(x){x}, ytrans, back = I),
    "2" = findbe2(bfct, efct, "Anke"),
    "3" = findbe3(),
    "4" = findbe2(bfct, efct, "Normolle"))
    
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        ## Finding initial values for c and d parameters
        cdVal <- findcd(x, y)
#        if (useFixed) {}  # not implemented at the moment
    
        ## Finding initial values for b and e parameters    
        beVal <- findbe(x, y, cdVal[1], cdVal[2])            
        
        ## Finding initial value for f parameter
        fVal <- 1  # better choice than 1 may be possible!   
    
        return(c(beVal[1], cdVal, beVal[2], fVal)[is.na(fixed)])
    }
}
