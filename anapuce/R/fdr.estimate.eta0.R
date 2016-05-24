fdr.estimate.eta0 <-
function(p,
    method=c("conservative", "adaptive", "bootstrap", "smoother"),
    lambda=seq(0,0.95,0.05) )
{
    method <- match.arg(method)
    
    # conservative method is default to force people
    # to think about their choice ...
    
    
    ########
    
    if (method == "conservative") # Benjamini and Hochberg (1995)
    {
        return(1.0)
    }

    ########
   
    if (method == "adaptive") # Benjamini and Hochberg (2000)
    {
        m <- length(p)
    sortp <- sort(p)
    s <- sort(1 - sortp)/(1:m)
    
    m0raw <- m
        i <- m
        while(i > 1 && s[i] <= s[i - 1]) i <- i - 1
        if(i > 1)
            m0raw <- 1/s[i - 1]
        else m0raw <- 1/s[1]
        
    m0 <- min(floor(1 + m0raw), m)
        
    eta0 <- m0/m   

        return(eta0)   
    }

    ########
    
    # for the remaining methods we require a lambda vector
    if (length(lambda)<4)
        stop("At least 4 values in lambda tuning vector required")
    
    eta0 <- rep(0,length(lambda))
    for(i in 1:length(lambda))
    {
        eta0[i] <- mean(p >= lambda[i])/(1-lambda[i])
    }
        
    ########
   
    if(method == "bootstrap") # Storey (2002) JRSSB
    {
            m <- length(p)
        mineta0 <- min(eta0)
            mse <- rep(0,length(lambda))
            eta0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    eta0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (eta0.boot-mineta0)^2
            }
            eta0 <- min(eta0[mse==min(mse)])
            eta0 <- min(eta0,1)
        
        return (eta0)
    }    
 
    ########
 
    if(method == "smoother") # Storey and Tibshirani (2003) PNAS
    {
            seta0 <- smooth.spline(lambda,eta0,df=3)
            eta0 <- predict(seta0,x=max(lambda))$y
            eta0 <- min(eta0,1)
        
        return(eta0)
    }   
 
}

