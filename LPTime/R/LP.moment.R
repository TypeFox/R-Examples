LP.moment <- function(x,m = 15){
    Tx <- LPTrans(x,m)
    c(mean(x),cor(x,Tx))
    }
