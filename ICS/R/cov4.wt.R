cov4.wt <- function(x, wt = rep(1/nrow(x), nrow(x)), location = TRUE, method="ML", na.action = na.fail)
    {
    if (length(wt) != nrow(x)) 
            stop("length of 'wt' must equal the number of rows in 'x'")
    method <- match.arg(method,c("ML","unbiased"))
    X <- data.frame(wt=wt)
    X$x <- x
    X <- na.action(X)
    x <- as.matrix(X$x)
    wt<- X$wt
    
    n <- nrow(x)
    p <- ncol(x)
    
        
        if (any(wt < 0) || (sum(wt)) == 0) 
            stop("weights must be non-negative and not all zero")
        
    
    COV.WT <- cov.wt(x, wt=wt, center=location, method=method)
    maha.wt<-mahalanobis(x,COV.WT$center, COV.WT$cov)
    if (is.logical(location)) 
        {
        location <- if (location) 
            COV.WT$center
        else 0
        }
    else {
        if (length(location) != p) 
            stop("length of 'location' must equal the number of columns in 'x'")
         }
    x.cent.weighted <- sqrt(wt) * sqrt(maha.wt) * sweep(x, 2, location, check.margin = FALSE)
    cov4.wt <-  crossprod(x.cent.weighted)/ ((p+2)* sum(wt))
    cov4.wt
    }
