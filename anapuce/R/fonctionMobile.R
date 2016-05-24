fonctionMobile <-
function(x, FUN, lag) 
{ 
        FUN=match.fun(FUN) 
        y<-numeric(length(x))
        n <- length(x) - lag
        # case 1: moving average of order lag
        y0 <- numeric(n)
        for (i in 1:n) { 
                y0[i] <- FUN(x[i:(i+lag)]) 
        } 
        # case 2: first values at left
        y1 <- numeric(ceiling(lag/2))
        for (i in 1:ceiling(lag/2)) {
        y1[i] <- FUN(x[1:(2*i-1)])
        }
        # The last values are considered equals to the last value calculated by mm
        y2 <- numeric(floor(lag/2))
        nmax <- length(x)
        nmin <- nmax-floor(lag/2)
        for (i in (nmin+1):nmax) {
        y2[nmax-i+1] <- y0[length(y0)]
        }
        y<-c(y1,y0,rev(y2))
        return(y) 
}

