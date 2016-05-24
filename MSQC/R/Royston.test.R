Royston.test <-
function(data) 
    {
    ## This function finds the equivalent degrees of freedom and calculates
    ## the Royston (1992) test statistic
    ## data is the data set being analyzed
    ## an n x p matrix (n = sample size, p = number of coordinates)

    ## Returns the test statistic and the approximated p-value


    "find.coeff.new"<-function(n)
    {
    ## This function finds the coefficients used to obtain the values after applying Royston's normalizing transformation.
        ## n is the number of observations
        ## p is the number of variables
        if(n <= 3) {
            stop(" n is too small!")
        }
        else if(n >= 4 && n <=11) {
            x = n
        }
        else if (n >= 12 && n <=2000){
            x = log(n)
        }
        else{
            stop("n is too large!")
        }


        if(n >= 4 && n <=11) {
            gamma = -2.273 + 0.459 * (x)
        }
        else if(n >= 12 && n <= 2000) {
            gamma = 0

        }
        else {
            stop("n is not in the proper size range")
        }
        ret <- list(gamma = gamma)
        if(n >= 4 && n <= 11) {
            mu = 0.5440 - 0.39978 * (x) + 0.025054 * (x^2) - 0.0006714 * (x^3)
        }
        else if(n >= 12 && n <= 2000) {
            mu = -1.5861 - 0.31082 * (x) - 0.083751 * (x^2) + 0.0038915 * (x^3)
        }
        else {
            stop(" n is not in the proper size range")
        }
        ret$mu <- mu
        if(n >= 4 && n <= 11) {
            sigma = exp(1.3822 - 0.77857 * (x) + 0.062767 * (x^2) - 0.0020322 * (x^3))
        }
        else if(n >= 12 && n <= 2000) {
            sigma = exp(-0.4803 -0.082676 * (x) + 0.0030302 * (x^2))
        }
        else {
            stop("n is not in the proper size range")
        }
        ret$sigma <- sigma
        return(ret)
    }

    "our.test.shapiro"<-function(data, n, p)
    {
    ## This function finds the values of Zj after the Royston transformation is applied
        ## data is the data file being analyzed, it is in matrix form
        ## n is the number of observations
        ## is the number of variables

        W <- rep(0, p)

        for(i in 1:p) {
            W[i] <- shapiro.test(data[, i])$statistic
        
}
        return(W)
    }



    "r.trans.new"<-function( n, p, a, b)
    {
    ## This function normalized the Shapiro-Wilks statistic under Royston's Transformation
        ## data is the data set used, in matrix form
        ## n is the number of observations
        ## p is the number of variables

        Z <- rep(0, p)
        if (n >= 4 && n <= 11)
        {
        for(i in 1:p) {
            Z[i] <- ((-log(a$gamma-(log(1-b[i])))-a$mu)/a$sigma)
        }
        }
        else if (n >= 12 && n <= 2000)
        {
            for (i in 1:p){
                Z[i] <- (((log(1-b[i])) + a$gamma - a$mu)/ a$sigma)
            }

        }
        else {
            stop ("WARNING : n is not in the proper range")

        }
        return(Z)
    }


    "Royston"<-function( n, p, z)
    {
## This function finds the statistics to be compared for multivariate normality
        ## data is the data set to be analysed
        ## n is the number of observations
        ## p is the number of variables

        G <- rep(0, p)
        for(i in 1:p) {
            G[i] <- (qnorm((pnorm( - z[i]))/2))^2
        }
        return(G)
    }

    d <- dim(data)
    n <- d[1]
    p <- d[2]
    a <- find.coeff.new(n)
    b <- our.test.shapiro(data, n, p)
    z <- r.trans.new( n, p, a, b)
    r <- Royston( n, p, z)
    u <- 0.715
    v <- 0.21364+0.015124*((log(n))^2)-0.0018034*((log(n))^3)
    la <- 5
    corr.data <- cor(data)
    new.corr <- ((corr.data^la) * (1 - (u  * (1 - corr.data)^
        u)/v))
    total <- sum(new.corr) - p
    avg.corr <- total/(p^2 - p)
    est.e <- p/(1 + (p - 1) * avg.corr)
    ans <- (est.e * (sum(r)))/p
    ## equivalent degrees of freedom
    #ans <- list(e = est.e)
    #ans <- list()
  test.statistic <- (est.e * (sum(r)))/p
    p.value <- 1 - (pchisq(test.statistic, df = est.e))
    #ans <- (est.e * (sum(r)))/p
    return(c("test.statistic"=test.statistic,"p.value"=p.value))
}
