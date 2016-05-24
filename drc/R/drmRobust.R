"drmRobust" <- function(robust, fctCall, lenData, lenPar)
{ 

    ## Finding robust scale estimate for trimmed and winsorised means and tukey
    if (robust%in%c("trimmed", "tukey", "winsor"))
    {
        call1 <- fctCall
        call1$"robust" <- "lts"

        scaleEst <- mad(residuals(eval(call1, parent.frame())), 0)
    }
    

    ## Defining distance functions
    quadratic <- function(x) {x*x}
    
    lms <- function(x) {median(x*x)}

    noRes <- floor((lenData+lenPar+1)/2)
    lts <- function(x) {sum(((x[order(x)])[1:noRes])^2)}
 
    metricTrim <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        c <- 1.345
        retVal <- x*x
        
        indexVec <- abs(x)>c
#        print(x)
        sumVec <- sum(indexVec)
#        print(sumVec)
        
        if (sumVec>0) {retVal[indexVec] <- rep(c*c, sumVec)}
        retVal
    }

    metricWinsor <- function(x) 
    {       
        if (all(is.na(x))) {return(x)}
    
        c <- 1.345

#        print(mad(x,0))
#        scaleEst <- (median(abs(x))/0.6745)  # overrules general scale estimate

#        scaleEst <- 9.03  # 9.055 for phones data set

        x <- x/scaleEst

        retVal <- x*x

        indexVec <- abs(x)>=c
        sumVec <- sum(indexVec)
        
        if (sumVec>0) {retVal[indexVec] <- (c*(2*abs(x)-c))[indexVec]}

#        retVal[abs(x)>c] <- (c*(2*abs(x)-c))[abs(x)>c] 
        retVal
    }

    tukeyBiweight <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        R <- 4.685
        retVal <- (x^6)/(R^4) - 3*(x^4)/(R^2) + 3*x*x
        
        indexVec <- abs(x)>R
        sumVec <- sum(indexVec)
        
        if (sumVec>0) {retVal[indexVec] <- rep(R*R, sumVec)}
        retVal
    }


    ## Assigning objective function
    robustFct <- switch(robust, mean=quadratic, median=abs, trimmed=metricTrim, tukey=tukeyBiweight, winsor=metricWinsor, lms = lms, lts = lts)
    
    return(robustFct)    
}
