pfNonlinBS <- function(data, particles=500, plot=FALSE) {
    if (missing(data)) {
         warning("data argument contained no data, using data simulated from the model.")
         data <- simNonlin(len=50)$data
    }
    res <- .Call("pfNonlinBS", data, particles, package="RcppSMC")

    time <- 1:length(data);
    if (plot) {
        with(res, plot(time,mean,type='l',
                       main='Filtering Mean and +/- 1,2 standard deviation intervals',
                       xlab='time',ylab='estimate',
                       xlim = c(0,length(data)),
                       ylim = c(min(mean-2.1*sd),max(mean+2.1*sd))))
        with(res,polygon(c(time,seq(length(data),1,-1)),
                         c(mean-2*sd,(mean+2*sd)[seq(length(data),1,-1)]),
                         col='palegreen1',border=NA))
        with(res,polygon(c(time,seq(length(data),1,-1)),
                         c(mean-1*sd,(mean+1*sd)[seq(length(data),1,-1)]),
                         col='palegreen3',border=NA))
        with(res,lines(time,mean, lwd=2, col='darkblue'))
    }

    invisible(res)
}



