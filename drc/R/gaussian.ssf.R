"gaussian.ssf" <- function(method = c("1", "2", "3", "4"), fixed, logg = FALSE, useFixed = FALSE)
{
    method <- match.arg(method)
        
    function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]

        ## Finding initial values for c and d parameters
        cdVal <- findcd(x, y)
        if (useFixed) {}  # not implemented at the moment        
    
        ## Finding initial values for b, e, and f parameters 
        if (logg) 
        {
            bVal <- 0.75 * sd(log(x[y > quantile(y, .75)]))
        } else {
            bVal <- 0.75 * sd(x[y > quantile(y, .75)])
        }   
        befVal <- c(bVal, x[which.max(y)], 1)
#        befVal <- c(sd(x), mean(x), 1)

        return(c(befVal[1], cdVal, befVal[2:3])[is.na(fixed)])
    }
}