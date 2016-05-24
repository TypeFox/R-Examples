############################################################################
### ------------------------------------------------
### In this comment substitute 'xxx' by
###  'var', 'median', 'IQR', 'mad',  respectively
### ------------------------------------------------
### We intentionally mask function 'xxx' from stats in order to add a formal
### argument "...".
### functionality of 'stats::xxx' is completely retained, however
### for help to the original 'stats::xxx' function write
###       'help("xxx", package="stats")'
### for code to the original 'stats::xxx' function write
###       'stats::xxx'
############################################################################


## ---- registering original functions as Methods: ----

setMethod("var", "ANY", function(x , ...)
       {dots <- list(...)
        if(hasArg(y)) y <- dots$"y"
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        if(!hasArg(use))
             use <- ifelse (na.rm, "complete.obs","all.obs")
        else use <- dots$"use"
        if(hasArg(y))
           stats::var(x = x, y = y, na.rm = na.rm, use)
        else
           stats::var(x = x, y = NULL, na.rm = na.rm, use)
        })

setMethod("median","ANY",function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        stats::median(x = x, na.rm = na.rm)}
)

setMethod("IQR","ANY",function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        stats::IQR(x = x, na.rm = na.rm)}
)

setMethod("mad","ANY",function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        low <- ifelse(hasArg(low), dots$"low", FALSE)
        high <- ifelse(hasArg(high), dots$"high", FALSE)
        center <- ifelse(hasArg(center), dots$"center", median(x, na.rm = na.rm))
        constant <- ifelse(hasArg(constant), dots$"constant", 1.4826)
        stats::mad(x = x, center = center, constant = constant, na.rm = na.rm, low = low, high = high)}
)


############################################################################
### acknowledgement:
###     methods 'skewness' and 'kurtosis' for particular methods
###     have been provided by G. Jay Kerns, gkerns@ysu.edu
###
### replaced e1071 version by sample and bias free (under normal distribution) 
### version for skewness and kurtosis (MK, 13 Nov. 2008)
############################################################################

setMethod("skewness","ANY",function(x , ...)
       {dots <- list(...)
        na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        sample.version <- ifelse(hasArg(sample.version), dots$"sample.version", FALSE)
        if (na.rm) x <- x[!is.na(x)]
        M <- mean(x)
        m3 <- mean((x-M)^3)
        m2 <- mean((x-M)^2)
        bias.cor <- ifelse(sample.version, 1, sqrt(n*(n-1))/(n-2))
        bias.cor*m3/m2^(3/2)
        })

setMethod("kurtosis","ANY",function(x , ...)
       {dots <- list(...)
        na.rm     <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
        sample.version <- ifelse(hasArg(sample.version), dots$"sample.version", FALSE)
        if (na.rm) x <- x[!is.na(x)]
        M <- mean(x)
        m4 <- mean((x-M)^4)
        m2 <- mean((x-M)^2)
        n <- length(x)
        g2 <- m4/m2^2 - 3
        if(sample.version){
            return(g2)
        }else{
            ## bias free for normal distributed samples
            return((n-1)/((n-2)*(n-3))*((n+1)*g2 + 6))
        }})
