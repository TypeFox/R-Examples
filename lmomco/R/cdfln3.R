"cdfln3" <-
function(x,para) {
    if(! are.parln3.valid(para)) return()
    XI <- para$para[1]
    U  <- para$para[2]
    A  <- para$para[3]

    f <- sapply(1:length(x), function(i) {
                               if(x[i] <= 0) return(NA)
                               Y <- (log(x[i]-XI) - U)/A
                               return(pnorm(Y)) })
    names(f) <- NULL
    return(f)
}

