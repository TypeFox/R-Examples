hzTest <-
function (data, cov = TRUE, qqplot = FALSE)
    {

        if (dim(data)[2] < 2 || is.null(dim(data))) {stop("number of variables must be equal or greater than 2")}
        if (!is.data.frame(data) && !is.matrix(data)) stop('Input must be one of classes \"data frame\" or \"matrix\"')
        
        dataframe=as.data.frame(data)
        dname <- deparse(substitute(data))
        data <- as.matrix(data)
        n <- dim(data)[1]
        p <- dim(data)[2]
        data.org = data
        
        if (cov){ 
            S <- ((n-1)/n)*cov(data)
        }
        else    {  
            S <- cov(data)
        }
        
        dif <- scale(data, scale = FALSE)
        
        
        Dj <- diag(dif%*%solve(S)%*%t(dif))  #squared-Mahalanobis' distances
        
        Y <- data%*%solve(S)%*%t(data)
        
        
        Djk <- - 2*t(Y) + matrix(diag(t(Y)))%*%matrix(c(rep(1,n)),1,n) + matrix(c(rep(1,n)),n,1)%*%diag(t(Y))
        
        b <- 1/(sqrt(2))*((2*p + 1)/4)^(1/(p + 4))*(n^(1/(p + 4))) #smoothing
{                                                                 #parameter
            if (qr(S)$rank == p){    
                HZ = n * (1/(n^2) * sum(sum(exp( - (b^2)/2 * Djk))) - 2 *
                              ((1 + (b^2))^( - p/2)) * (1/n) * (sum(exp( - ((b^2)/(2 *
                                                                                       (1 + (b^2)))) * Dj))) + ((1 + (2 * (b^2)))^( - p/2)))
            }
            else {
                HZ = n*4
            }  
            
        }
        wb <- (1 + b^2)*(1 + 3*b^2)
        
        a <- 1 + 2*b^2
        
        mu <- 1 - a^(- p/2)*(1 + p*b^2/a + (p*(p + 2)*(b^4))/(2*a^2)) #HZ mean
        
        si2 <- 2*(1 + 4*b^2)^(- p/2) + 2*a^( - p)*(1 + (2*p*b^4)/a^2 + (3*p*
                                                                            (p + 2)*b^8)/(4*a^4)) - 4*wb^( - p/2)*(1 + (3*p*b^4)/(2*wb) + (p*
                                                                                                                                               (p + 2)*b^8)/(2*wb^2)) #HZ variance
        
        pmu <- log(sqrt(mu^4/(si2 + mu^2))) #lognormal HZ mean
        psi <- sqrt(log((si2 + mu^2)/mu^2)) #lognormal HZ variance
        
        P <- 1 - plnorm(HZ,pmu,psi) #P-value associated to the HZ statistic
        
        
        
        if (qqplot){    
            d <- Dj    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot",
                 xlab = "Squared Mahalanobis Distance",ylab="Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }
        
        
        result <- new("hz", HZ = HZ, p.value = P, dname = dname, dataframe = dataframe)
        
        result
        
    }
