##############################
## Functions to plot things ##
##############################

plot.true.functions <- function(x, y, z,
                                f1, f2){
    plot(x, y,
         ylab='f(x)', xlab='x', main='True functions')
    lines(x[z==1], y[z==1], pch=21, col=1, type='p')
    lines(x[z==2], y[z==2], pch=21, col=2, type='p')
    lines(x, f1)
    lines(x, f2, col=2)

    legend('bottomleft',
           c('f_1 true', 'f_2 true'),
           col=c(1, 2),
           lty=c(1, 1))
}


plot.initial.functions <- function(x, y, z,
                                   f1.true, f2.true,
                                   f1.initial, f2.initial){
    plot(x, y,
         ylab='f(x)', xlab='x', main='Initial functions')
    lines(x[z==1], y[z==1], pch=21, col=1, type='p')
    lines(x[z==2], y[z==2], pch=21, col=2, type='p')
    lines(x, f1.true)
    lines(x, f2.true, col=2)
    lines(x, f1.initial, col=1, lty=2)
    lines(x, f2.initial, col=2, lty=2)
    
    legend('bottomleft',
           c('true f_2', 'initial f_2', 'true f_1', 'initial f_1'),
           col=c(2, 2, 1, 1),
           lty=c(1, 2, 1, 2))
}


plot.estimated.functions <- function(x, y, z, p1k,
                                     f1.true, f2.true,
                                     f1hat, f2hat){
    zhat <- rep(1, length(x))
    zhat[which(p1k < 0.5)] = 2
    
    plot(x, y, pch=20,
         ylab='f_i(x)', xlab='x')
    lines(x[z==1], y[z==1], col=1, pch=20, type='p', ylab='f(x)', xlab='x')
    lines(x[z==2], y[z==2], col=2, pch=20, type='p')
    lines(x[which(z!=zhat)], y[z!=zhat], col=3, pch=22, type='p')
    lines(x, f1.true)
    lines(x, f2.true, col=2)
    lines(x, f2hat, lty=2, col=2)
    lines(x, f1hat, lty=2, col=1)
    
    legend('bottomleft',
           c('true f_2', 'f_2hat', 'true f_1', 'f_1hat'),
           col=c(2, 2, 1, 1),
           lty=c(1, 2, 1, 2))
}
