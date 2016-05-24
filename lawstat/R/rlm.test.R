`rlm.test` <-
function (x, crit.values = c("chisq.approximation", "empirical"), N = 0) 
{

### set the default option to be chi-square approximation ###
    
    crit.values = match.arg(crit.values)

### stop the function if it is not a vector ###

    if (NCOL(x) > 1) {
        stop("x is not a vector or univariate time series")
    }

### stop the function if the data has missing values ###

    if (any(is.na(x))) {
        stop("NAs in x")
    }

### stop the function if the number of Monte Carlo simulations is not provided ###

    if ((crit.values == "empirical") & (N == 0)) {
        stop("number of Monte Carlo simulations N should be provided for the empirical critical values")
    }

### initialize variables ###

    DNAME <- deparse(substitute(x))
    n <- length(x)
    m1 <- sum(x)/n
    m3 <- sum((x - m1)^3)/n
    m4 <- sum((x - m1)^4)/n

### construct the test statistic ###

    J <- sqrt(2)*mean(abs(x - median(x)))
    ek <- 6
    b1 <- (m3/(J^3))^2
    b2 <- (m4/(J^4)-ek)^2
    vk <- 1200/n
    vs <- 60/n
    statistic <- b1/vs + b2/vk
    METHOD <- "Robust L1 moment-based goodness-of-fit test"

### calculate empirical critical values ###
    
    if (crit.values == "empirical") 
    {

    METHOD <- paste(METHOD, "using empirical critical values with N =", N)

### create a vector "jb" to store statistics ###

    jb <- double(N)

### generate random Laplace variables to calculate critical values ###

        for (k in 1:N) 
        {
        e <- rexp(n)-rexp(n)
                m1 <- sum(e)/n
    m3 <- sum((e - m1)^3)/n
    m4 <- sum((e - m1)^4)/n
    J <- sqrt(2)*mean(abs(e - median(e)))
    ek <- 6
    b1 <- (m3/(J^3))^2
    b2 <- (m4/(J^4)-ek)^2
                vk <- 1200/n
                vs <- 60/n
                jb[k] <- b1/vs + b2/vk
            }

### sort the generated statistics ###

            y <- sort(jb)

### set the p-value to zero if the statistic is greater than maximum of generated statistics ###

            if (statistic >= max(y)) 
    {
                p.value = 0
            }

### set the p-value to one if the statistic is smaller than minimum of generated statistics ###

            else if (statistic <= min(y)) 
    {
                p.value = 1
            }

### calculate the p-value in the case the statistic is between min and max of generated statistics ###

            else 
    {
                an <- which(y == max(y[I(y < statistic)]))
                bn <- which(y == min(y[I(y >= statistic)]))
                a <- max(y[I(y < statistic)])
                b <- min(y[I(y >= statistic)])
                pa <- (an - 1)/(N - 1)
                pb <- (bn - 1)/(N - 1)
                alpha <- (statistic - a)/(b - a)
                p.value = 1 - alpha * pb - (1 - alpha) * pa
            }
        
    }

### calculate the p-value using a Chi-squared approximation ###

    else 
    {
METHOD <- paste(METHOD,"using a Chi-squared approximated critical values")
p.value <- pchisq(statistic, df = 2, lower.tail = FALSE)
    }

### display output ###
    
    STATISTIC = statistic
    names(STATISTIC) <- paste("Chi-squared statistic")
    PARAMETER <- 2
    names(PARAMETER) <- "df"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = p.value, method = METHOD, data.name = DNAME), 
        class = "htest")
}

