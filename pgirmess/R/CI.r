"CI"<-function(x, y, totrials=FALSE){
    if (totrials) totpop<-y else totpop<-x+y
    bashau <- cbind(rep(0, length(x)), rep(0, length(x)))
    for(i in 1:length(x)) {
        ic <- prop.test(x[i], totpop[i])$conf.int
        bashau[i, 1] <- ic[1]
        bashau[i, 2] <- ic[2]
    }
    bashau
}
