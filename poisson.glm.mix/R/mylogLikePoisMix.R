mylogLikePoisMix <-
function(y, mean, pi){
    if (is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
        stop(paste(sQuote("y"), "must be a matrix"))
    if (min(y) < 0 | sum(round(y)) != sum(y)) 
        stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
    if (length(mean) != length(pi)) 
        stop(paste(sQuote("mean"), "must be a list of the same length as", 
            sQuote("pi")))
    g <- length(pi)
    n <- dim(y)[1]
    cols <- dim(y)[2]
    nas<-0
    y <- matrix(y, nrow = n, ncol = cols)
    logLike <- rep(0, n)
index<-1:g
    epsilon<-exp(-720)
thresh<- -745
nn<-0
logpi<-log(pi)
#for (i in 1:n) {
ef <- matrix(logpi,nrow=n,ncol=g,byrow=T)
for (k in 1:g){
ef[,k] <- ef[,k] + rowSums(dpois(y,mean[[k]],log = T))
}
efmax<-apply(ef,1,max)
ef <- ef - efmax
logLike <- efmax + log(rowSums(exp(ef)))
#}
#if (nn>0) print(paste("warning: for ", nn, " observations, an upper bound of the loglikelihood was computed"))
#print(paste("warning:", nas, " NAs"))
        
    return(list(ll = sum(logLike, na.rm = TRUE)))
}
