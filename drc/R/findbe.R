## Defining functions for finding initial values of the b and e parameter

## Finding b and e based on linear regression after logit transformation (classical approach)
findbe1 <- function(doseTr, respTr, sgnb = 1, back = exp)
{
    function(x, y, cVal, dVal)
    {
#        lmFit <- lm(log((dVal - y)/(y - cVal)) ~ log(x), subset = x > 0)

#respTr <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
#doseTr <- function(x) {rVec <- log(x); rVec[!x>0] <- NA; rVec}
        lmFit <- lm(respTr(y, cVal, dVal) ~ doseTr(x))
        coefVec <- coef(lmFit)
        bVal <- sgnb * coefVec[2]        
        eVal <- back(-coefVec[1] / (sgnb * bVal))

        return(as.vector(c(bVal, eVal)))
    }
}  

## Anke's procedure
findbe2 <- function(bfct, efct, method, sgnb = 1)
{
    ## Helper functions used below
    # bfct <- function(x, y, cVal, dVal, eVal) {log((dVal - y) / (y - cVal)) / log(x / eVal)}
    bFct <- function(x, y, cVal, dVal, eVal)
    {
        median(bfct(x, y, cVal, dVal, eVal), na.rm = TRUE)
    }
    # efct <- function(x, y, bVal, cVal, dVal) {x * (((dVal - y) / (y - cVal))^(-1 / bVal))}
    eFct <- function(x, y, bVal, cVal, dVal)
    {
        median(efct(x, y, bVal, cVal, dVal), na.rm = TRUE)
    }
    
    switch(method,
    "Anke" = 
    function(x, y, cVal, dVal)
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
#        print(maxDose)
        
        ## Smallest dose with all responses below
        belowVec <- x[y < midResp]
        uniBelow <- unique(belowVec)
        if (length(belowVec) < sum(x %in% uniBelow))
        {
            uniBelow <- tail(uniBelow, -1)
        }
        minDose <- min(x[x %in% uniBelow])
#        print(minDose)

        subsetInd <- (x > maxDose) & (x < minDose)
#        print(subsetInd)
        eVal <- mean((1 / (1 + (abs(y[subsetInd] - midResp)/15)^2)) * x[subsetInd])
        if (is.nan(eVal))  # in case subsetInd are all FALSE
        {
            eVal <- (minDose + maxDose) / 2
        }
        sort1 <- sort(unique(x))[2]
        if (eVal < sort1)
        {
            eVal <- sort1
        }
        sort2 <- sort(unique(x), decreasing = TRUE)[2]
        if (eVal > sort2)
        {
            eVal <- sort2
        }
#        print(eVal)

        ## Finding initial value for b
#        bVal <- median(log((dVal - y) / (y - cVal)) / log(x / eVal), na.rm = TRUE)
        bVal <- bFct(x, y, cVal, dVal, eVal)
#        print(bVal)
        
        ## Checking sign of b and possibly take action if it's wrong
        regSlope <- as.vector(coef(lm(y ~ x)))[2]
        if ((!is.na(bVal)) && ((sgnb * regSlope / bVal) > 0))
        {
#             bVal <- -bVal
#             eVal <- median(x * (((dVal - y) / (y - cVal))^(-1/bVal)), na.rm = TRUE)
             eVal <- eFct(x, y, -bVal, cVal, dVal)
             bVal <- bFct(x, y, cVal, dVal, eVal)   
        }
        if (is.na(bVal))
        {
            bVal <- (sgnb) * (-regSlope)
        }
        print(c(bVal, eVal))
        return(c(bVal, eVal))
    },
    "Normolle" = 
    function(x, y, cVal, dVal)
    {
        bVal <- bFct(x, y, cVal, dVal, mean(range(x)))
        eVal <- eFct(x, y, bVal, cVal, dVal)

        return(c(bVal, eVal)) 
    })    
}
    
## Finding b and e based on stepwise increments
findbe3 <- function(sgnb = 1)
{
    function(x, y, cVal, dVal)
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
        bVal <- (sgnb) * (-sign(uniy[j] - uniy[j-1]))  # -(uniy[j] - uniy[j-1]) / (unix[j] - unix[j-1])

        return(as.vector(c(bVal, eVal)))  
    }
}
    
    
## Normolle's procedure
#findbe4 <- function(bfct, efct)
#{
#    ## Helper functions used below
#    # bfct <- function(x, y, cVal, dVal, eVal) {log((dVal - y) / (y - cVal)) / log(x / eVal)}
#    bFct <- function(x, y, cVal, dVal, eVal)
#    {
#        median(bfct(x, y, cVal, dVal, eVal), na.rm = TRUE)
#    }
#    # efct <- function(x, y, bVal, cVal, dVal) {x * (((dVal - y) / (y - cVal))^(-1 / bVal))}
#    eFct <- function(x, y, bVal, cVal, dVal)
#    {
#        median(efct(x, y, bVal, cVal, dVal), na.rm = TRUE)
#    }
#
#    function(x, y, cVal, dVal)
#    {
##        initeVal <- mean(range(x))
#        bVal <- bFct(x, y, cVal, dVal, mean(range(x)))
#        eVal <- eFct(x, y, bVal, cVal, dVal)
#
#        return(c(bVal, eVal)) 
#    }    
#}
