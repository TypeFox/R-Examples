`Burstsr` <-
function(x){

y<-Bursts(x)

y <- c(y, rev(y))
        n <- length(y)
    y[seq(from=1, to=n, by=2)]


}

