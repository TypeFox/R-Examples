epi.stratasize <- function (strata.n, strata.mean, strata.var, strata.Py, epsilon.r, method = "mean", conf.level = 0.95) 
{
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    if (method == "total") {
        N <- sum(strata.n)
        mean <- sum(strata.n * strata.mean) / N
        sigma.bx <- sum(strata.n * (strata.mean - mean)^2) / N
        sigma.wx <- sum(strata.n * strata.var) / N
        sigma.x <- sigma.bx + sigma.wx
        V <- sigma.x/mean^2
        gamma <- sigma.bx/sigma.wx
        
        # Equation 6.25 Levy and Lemeshow. Example on p 177 gives 9 for z^2 which equates to an alpha of a bit less than 0.01. 
        # Emailed Stan Lemeshow re this on 9 Feb 2008 - he confirms this is true.
        total.sample <- round((((z^2 * N)/(1 + gamma)) * V) / (((z^2 * V) / (1 + gamma)) + N * (epsilon.r^2)), digits = 0)
        strata.sample <- round(strata.n * (total.sample/N), digits = 0)
        
        result.01 <- c(strata.sample)
        result.02 <- c(total.sample)
        result.03 <- cbind(mean = mean, sigma.bx = sigma.bx, 
            sigma.wx = sigma.wx, sigma.x = sigma.x, rel.var = V, 
            gamma = gamma)
        
        rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    }
    
    if (method == "mean") {
        N <- sum(strata.n)
        mean <- sum(strata.n * strata.mean) / N
        sigma.bx <- sum(strata.n * (strata.mean - mean)^2) / N
        sigma.wx <- sum(strata.n * strata.var) / N
        sigma.x <- sigma.bx + sigma.wx
        V <- sigma.x/mean^2
        gamma <- sigma.bx/sigma.wx
        
        # Equation 6.25 Levy and Lemeshow. Example on p 177 gives 9 for z^2. Suspect this is an error. I use 1.96^2 =~ 4
        total.sample <- round((((z^2 * N)/(1 + gamma)) * V) / (((z^2 * V) / (1 + gamma)) + N * (epsilon.r^2)), digits = 0)
        strata.sample <- round(strata.n * (total.sample/N), digits = 0)
        
        result.01 <- c(strata.sample)
        result.02 <- c(total.sample)
        result.03 <- cbind(mean = mean, sigma.bx = sigma.bx, 
            sigma.wx = sigma.wx, sigma.x = sigma.x, rel.var = V, 
            gamma = gamma)
        
        rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    }
    
    if (method == "proportion") {
        # Where method == "proportion" the estimated proportions for each strata are entered into the vector strata.Py:
        N <- sum(strata.n)
        mean <- sum(strata.n * strata.Py) / N
        # The vector strata.var is ignored (variance of proportion calculated as follows):
        strata.var = (strata.Py * (1 - strata.Py))
        phi <- (strata.n * sqrt(strata.var))/sum(strata.n * sqrt(strata.var))
        sigma.bx <- sum((strata.n^2 * strata.var)/((phi) * (mean^2)))
        sigma.bxd <- sum((strata.n * strata.var)/mean^2)
        
        # Equation 6.23 Levy and Lemeshow. Note the similarity between 6.23 and 6.22:
        total.sample <- round(((z^2/N^2) * sigma.bx)/((epsilon.r^2) + ((z^2/N^2) * sigma.bxd)), digits = 0)
        strata.sample <- round(strata.n * (total.sample/N), digits = 0)
        
        result.01 <- c(strata.sample)
        result.02 <- c(total.sample)
        result.03 <- cbind(mean = mean, sigma.bx = sigma.bx, sigma.bxd = sigma.bxd, phi = phi)
        
        rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    }
    
    if (method == "pps") {
        N <- sum(strata.n)
        mean <- sum(strata.n * strata.mean)/N
        strata.var = (strata.mean * (1 - strata.mean))
        sigma.bx <- sum((strata.n * strata.var)/mean^2)
        
        total.sample <- round(((z^2/N) * sigma.bx)/(epsilon.r^2 + ((z^2/N^2) * sigma.bx)), digits = 0)
        strata.sample <- round(strata.n * (total.sample/N), digits = 0)
        
        result.01 <- c(strata.sample)
        result.02 <- c(total.sample)
        result.03 <- cbind(mean = mean, sigma.bx = sigma.bx)
        
        rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    }
    return(rval)
}
