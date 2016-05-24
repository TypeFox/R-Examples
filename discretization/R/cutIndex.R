cutIndex <-
function(x,y){
       n <- length(y)
       initEnt <- 9999 
       entropy <- initEnt; ci <- NULL;
       for (i in 1:(n-1)){
            if(x[i+1]!=x[i]) {
                ct <- (x[i]+x[i+1])/2
                wx <- which(x<=ct)
                wn <- length(wx)/n
                e1 <- wn*ent(y[wx])
                e2 <- (1-wn)*ent(y[-wx])
                val <- e1+e2
                if(val<entropy) {
                    entropy <- val
                    ci <- i
                }
            }
       }
       if(is.null(ci)) return(NULL) 
       return (c(ci, entropy))
}
