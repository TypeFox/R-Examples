HZ.test <-
function(data) 
    {
    ## This function calculates the Henze-Zirkler test
    ## statistic for assessing  multivariate normality of a data set to be analyzed.
    ## data is the data set to be analyzed.
    ## r is the desired alpha level of significance
    ## returns the test statistic and the approximated p-value
    data <- as.matrix(data);
    "loop4"<- function(data)
    {
    #data<-as.matrix(data); #YO se lo agregue
    n <- dim(data)[1]
    x <- data
    s.inv <- solve(var(x))
    X <- X.2 <- data %*% s.inv %*% t(data)
    diag(X) <- 0
    for(i in 1:(n-1)) {
        for(j in (i+1):n) {
            X[j,i] <- X[i, j] <- X.2[i,i] - 2 * X.2[i,j] + X.2[j,j]
        }
    }
    return(X)
    }

    "Square"<-function(data)
    {
        S <- var(data)
        S.inv <- solve(S)
        x <- data
        x.diff <- scale(x, scale = F)
        square <- diag(x.diff %*% S.inv %*% t(x.diff))
        return(square)
    }


    d <- dim(data)
    n <- d[1]
    p <- d[2]


    Beta <- 1/(sqrt(2)) * ((2 * p + 1)/4)^(1/(p + 4)) * (n^(1/(p + 4)))

    a <- var(data)
    m <- dim(a)[1]
    ranka <- qr(a)$rank

    if( ranka == m){

    s <- loop4(data)
    square <- Square(data)

    T.beta <- n * (1/(n^2) * (sum(exp( - (Beta^2)/2 * s))) - 2 * ((1 + (
        Beta^2))^( - p/2)) * (1/n) * (sum(exp( - ((Beta^2)/(2 * (1 +
        (Beta^2)))) * square))) + ((1 + (2 * (Beta^2)))^( - p/2)))
    }
    else
    {
        T.beta <- n * 4
    }
    ans <- list(T = T.beta)
    mu <- 1 - ((1 + 2 * (Beta^2))^( - p/2)) * (1 + (p * (Beta^2))/(1 + 2 *
        (Beta^2)) + (p * (p + 2) * (Beta^4))/(2 * (1 + 2 * (Beta^2))^
        2))

    w.beta <- (1 + (Beta^2)) * (1 + 3 * (Beta^2))
    sigma.sq <- 2 * (1 + 4 * (Beta^2))^( - p/2) + 2 * (1 + 2 * (Beta^2))^
        ( - p) * (1 + (2 * p * (Beta^4))/(1 + (2 * (Beta^2)))^2 + (
        3 * p * (p + 2) * (Beta^8))/(4 * (1 + 2 * (Beta^2))^4)) - 4 *
        (w.beta^( - p/2)) * (1 + (3 * p * (Beta^4))/(2 * w.beta) + (
        p * (p + 2) * (Beta^8))/(2 * (w.beta)^2))


    p.mu<-log(sqrt((mu^4)/(sigma.sq + (mu^2))))

    p.sigma<-sqrt(log((sigma.sq + (mu^2))/(mu^2)))

  
    p.value<- 1 - (plnorm(T.beta,p.mu,p.sigma))
    test.statistic <- T.beta
    ans <- T.beta; #Yo se lo agregue
    return(c(p.value,test.statistic))
}
