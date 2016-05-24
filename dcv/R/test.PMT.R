test.PMT <- function(x, y){
    n <- length(x)
    ym1 <- c()
    ym2 <- c()
    pyy2 <- c()
    for(i in 1:n){
        ym1[i] <- x[i] - mean(x)
        ym2[i] <- y[i] - mean(y)
        pyy2[i] <- x[i] * y[i]
    }
    
    z <- 0
    l <- 0
    yy1 <- c()
    zz1 <- c()
    for(i in 1:n){
        if(pyy2[i] > 0){
            z <- z+1
            yy1[z] <- pyy2[i]
        }
        else{
            l <- l+1
            zz1[l] <- abs(pyy2[i])
        }
    }
    s1 <- 0
    for(i in 1:z){
        s1 <- s1+(yy1[i]-mean(yy1))^2
    }
    s1 <- s1/(z - 1)
    s2 <- 0
    for(i in 1:l){
        s2 <- s2+(zz1[i]-mean(zz1))^2
    }
    s2 <- s2/(l - 1)
    t  <- (mean(yy1)-mean(zz1))/sqrt(s1/z + s2/l)
    class(t) <- "PMT"
    return(t)
}

