chiSq <-
function(tb){## this works only for 2 x k tables
     tb <- tb+0.0001 
     e <- tb
     n <- sum(tb)
     p <- length(tb[,1]); q <- length(tb[1,])
     mi <- numeric(p); ni <- numeric(q)
     for(i in 1:q) ni[i] <- sum(tb[,i])
     for(i in 1:p) {
        mi[i] <- sum(tb[i,])
        e[i,] <- mi[i]*ni/n
     }     
     val <- sum((tb-e)^2/e)
     return(val)
}
