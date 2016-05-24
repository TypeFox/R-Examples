plot1.acf <-
function(x, xlab='Lag',type='h', ylab='ACF',ci.type=c("white","ma")[1],
main=paste('Series', x$series),ci.col='blue',  ...){
bound=switch(ci.type, 
               white=1.96/sqrt(x$n.used),
               ma={ wt <- sqrt(cumsum(c(1,2*x$acf^2)))
                               wt <- wt[-length(wt)]
                               max(1.96/sqrt(x$n.used)*wt)}
)
ylim <- c(min(c(min(x$acf), -bound)), max(c(max(x$acf), bound)))
plot(y=x$acf, x=x$lag,type=type, xlab=xlab, ylab=ylab, ylim=ylim,main=main,...)
abline(h=0)
switch(ci.type, 
               white={abline(h=1.96/sqrt(x$n.used), col=ci.col, lty=2)
                                 abline(h=-1.96/sqrt(x$n.used), col=ci.col, lty=2)},
               ma={ wt <- sqrt(cumsum(c(1,2*x$acf^2)))
                               wt <- wt[-length(wt)]
                               lines(y=1.96/sqrt(x$n.used)*wt, x=x$lag,col=ci.col,lty=2)
                               lines(y=-1.96/sqrt(x$n.used)*wt, x=x$lag,col=ci.col,lty=2)                                 
})
invisible() 
}
