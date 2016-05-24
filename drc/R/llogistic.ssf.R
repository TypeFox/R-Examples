"llogistic.ssf" <- function(method = c("1", "2", "3", "4"), fixed, useFixed = FALSE)
{
    method <- match.arg(method)
    
    ## Defining helper functions (used below)
    ytrans <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
    bfct <- function(x, y, cVal, dVal, eVal) {ytrans(y, cVal, dVal) / log(x / eVal)}
#    efct <- function(x, y, bVal, cVal, dVal) {x * (((dVal - y) / (y - cVal))^(-1 / bVal))}
    efct <- function(x, y, bVal, cVal, dVal) {x * exp(-ytrans(y, cVal, dVal)/bVal)}
    
    ## Assigning function for finding initial b and e parameter values    
    findbe <- switch(method,
    "1" = findbe1(function(x) {rVec <- log(x); rVec[!x>0 | !is.finite(x)] <- NA; rVec}, ytrans),
    "2" = findbe2(bfct, efct, "Anke"),
    "3" = findbe3(),
    "4" = findbe2(bfct, efct, "Normolle"))
    
    function(dframe)
    {
        ncoldf <- ncol(dframe)
        x <- dframe[, 1]        
#        x <- dframe[, -ncoldf]
        y <- dframe[, ncoldf]
    
#        x <- dframe[, 1]
#        y <- dframe[, 2]

        ## Finding initial values for c and d parameters
        cdVal <- findcd(x, y)
#        if (useFixed) {  # not implemented at the moment
#            cdVal <- c(ifelse(notFixed[2], cdVal[1], fixed[2] / respScaling),
#            ifelse(notFixed[3], cdVal[2], fixed[3] / respScaling))}
    
        ## Finding initial values for b and e parameters    
        beVal <- findbe(x, y, cdVal[1], cdVal[2])       
    
        ## Finding initial value for f parameter
        fVal <- 1
        # better choice than 1 may be possible! 
        # the f parameter, however, is very rarely a magnitude of 10 larger or smaller
        
        return(c(beVal[1], cdVal, beVal[2], fVal)[is.na(fixed)])
    }
}

if (FALSE)
{

"llogistic.ssf" <- function(method = c("1", "2", "3", "4"), fixed)
{
    method <- match.arg(method)
    numParm <- length(fixed)
    notFixed <- is.na(fixed)

    ## Version 1 (default)    
    ssFct <- switch(method,
    "1" = 
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        yRange <- range(y)
        lenyRange <- 0.001 * diff(yRange)
        cVal <- yRange[1] - lenyRange  # the c parameter        
        dVal <- yRange[2] + lenyRange  # the d parameter
        fVal <- 1  # better choice may be possible!        

        beVal <- find.be1(x, y, cVal, dVal)
        
        return(c(beVal[1], cVal, dVal, beVal[2], fVal)[notFixed])
    },
    "2" =
    function(dframe, doseScaling, respScaling)
    {
        x <- dframe[, 1] / doseScaling
        y <- dframe[, 2] / respScaling

        cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2] / respScaling)
        dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3] / respScaling)

        fVal <- 1  # need not be updated with value in 'fixed[5]'
        # better choice than 1 may be possible! 
        # the f parameter, however, is very rarely a magnitude of 10 larger or smaller

#        # Cutting away response values close to d
#        indexT1a <- x > 0
#        x2 <- x[indexT1a]
#        y2 <- y[indexT1a]
            
        beVal <- find.be1(x, y, cVal, dVal)
# These lines are not needed as the b and e parameters are not used in further calculations
#        bVal <- ifelse(notFixed[1], beVec[1], fixed[1])
#        eVal <- ifelse(notFixed[4], beVec[2], fixed[4] / doseScaling)

#        bVal <- beVec[1]
#        eVal <- beVec[2]                
#        return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
        
        return(c(beVal[1], cVal, dVal, beVal[2], fVal)[notFixed])
    },
    "3" =
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        cVal <- ifelse(notFixed[2], 0.99 * min(y), fixed[2])
        dVal <- ifelse(notFixed[3], 1.01 * max(y), fixed[3])
        fVal <- 1  # need not be updated with value in 'fixed[5]'
        
#        if ( length(unique(x)) == 1 ) {return((c(NA, NA, dVal, NA, NA))[notFixed])}  
#        # only estimate of upper limit if a single unique dose value 
# no longer needed
           
        beVec <- find.be2(x, y, cVal, dVal)
        bVal <- beVec[1]
        eVal <- beVec[2]
            
        return(as.vector(c(bVal, cVal, dVal, eVal, fVal)[notFixed]))
    })

    ## Finding b and e based on linear regression
    find.be1 <- function(x, y, cVal, dVal)
    {
#        logitTrans <- log((d - y)/(y - c))  
        lmFit <- lm(log((dVal - y)/(y - cVal)) ~ log(x), subset = x > 0)
        coefVec <- coef(lmFit)
        bVal <- coefVec[2]        
        eVal <- exp(-coefVec[1] / bVal)

        return(as.vector(c(bVal, eVal)))
    }  

    ## Finding b and e based on stepwise increments
    find.be2 <- function(x, y, cVal, dVal)
    {
        unix <- unique(x)
        uniy <- tapply(y, x, mean)
        lenx <- length(unix)
        
        j <- 2
        for (i in 2:lenx)
        {
            crit1 <- (uniy[i] > (cVal + dVal)/2) && (uniy[i-1] < (cVal + dVal)/2)
            crit2 <- (uniy[i] < (cVal + dVal)/2) && (uniy[i-1] > (cVal + dVal)/2)
            if (crit1 || crit2) break
            j <- j + 1
        }
        eVal <- (unix[j] + unix[j-1]) / 2
        bVal <- -sign(uniy[j] - uniy[j-1])  # -(uniy[j] - uniy[j-1]) / (unix[j] - unix[j-1])
        return(as.vector(c(bVal, eVal)))  
    }
    
#    ## Finding b and e based on linear regression
#    find.be2 <- function(x, y, c, d)
#    {
#        logitTrans <- log((d - y)/(y - c))  
#
#        lmFit <- lm(logitTrans ~ log(x))
##        eVal <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))
##        bVal <- coef(logitFit)[2]
#
#        coefVec <- coef(lmFit)
#        bVal <- coefVec[2]        
#        eVal <- exp(-coefVec[1]/bVal)    
#
#        return(as.vector(c(bVal, eVal)))
#    }
    
    ## Defining function for finding initial values of the b and e parameter
    
    ## Helper functions used below
    bFct <- function(x, y, cVal, dVal, eVal)
    {
        median(log((dVal - y) / (y - cVal)) / log(x / eVal), na.rm = TRUE)
    }

    eFct <- function(x, y, cVal, dVal, bVal)
    {
#        print((dVal - y) / (cVal - y))
        median(x * (((dVal - y) / (y - cVal))^(-1 / bVal)), na.rm = TRUE)
    }
    
    ## Anke's procedure
    find.be3 <- function(x, y, cVal, dVal)
    {
        ## Finding initial value for e
        midResp <- (dVal - cVal) / 2
        
        ## Largest dose with all responses above
        aboveVec <- x[y > midResp]
        uniAbove <- unique(aboveVec)
        
        if (length(aboveVec) < sum(x %in% uniAbove))
        {
            uniAbove <- head(uniAbove, -1)
        }
        maxDose <- max(x[x %in% uniAbove])
        print(maxDose)
        
        ## Smallest dose with all responses below
        belowVec <- x[y < midResp]
        uniBelow <- unique(belowVec)
        if (length(belowVec) < sum(x %in% uniBelow))
        {
            uniBelow <- tail(uniBelow, -1)
        }
        minDose <- min(x[x %in% uniBelow])
        print(minDose)

        subsetInd <- (x > maxDose) & (x < minDose)
        eVal <- mean((1 / (1 + (abs(y[subsetInd] - midResp)/15)^2)) * x[subsetInd])

        ## Finding initial value for b
#        bVal <- median(log((dVal - y) / (y - cVal)) / log(x / eVal), na.rm = TRUE)
        bVal <- bFct(x, y, cVal, dVal, eVal)
        
        ## Checking sign of b and possibly take action if it's wrong
        if ((coef(lm(y ~ x))[2] / bVal) > 0)
        {
#             bVal <- -bVal
#             eVal <- median(x * (((dVal - y) / (y - cVal))^(-1/bVal)), na.rm = TRUE)
             eVal <- eFct(x, y, cVal, dVal, -bVal)
             bVal <- bFct(x, y, cVal, dVal, eVal)   
        }
        
        return(c(bVal, eVal))
    }
    
    ## Normolle's procedure
    find.be4 <- function(x, y, cVal, dVal)
    {
        initeVal <- mean(range(x))
        bVal <- bFct(x, y, cVal, dVal, initeVal)
        eVal <- eFct(x, y, cVal, dVal, bVal)

#        bVal <- bFct(x, y, cVal, dVal, eVal) 
#        eVal <- eFct(x, y, cVal, dVal, bVal)

        return(c(bVal, eVal)) 
    }    

    ## Returning self starter function
    ssFct
}

LL.ssf <- llogistic.ssf
}

